from __future__ import annotations

import csv
import logging
from collections import defaultdict
from pathlib import Path

import networkx as nx

logger = logging.getLogger(__name__)

STRONG_EVIDENCE_TYPES = {
    "physical_overlap",
    "segment_membership",
    "co_abundance_guild",
}


class BinningMiner:
    def __init__(
        self,
        graph: nx.MultiGraph,
        strong_evidence_types: set[str] | None = None,
        min_rescue_support: int = 1,
        dense_membership: bool = False,
    ):
        self.graph = graph
        self.strong_evidence_types = strong_evidence_types or set(STRONG_EVIDENCE_TYPES)
        self.min_rescue_support = min_rescue_support
        self.dense_membership = dense_membership

    def run(self, binning_tsv: str | Path) -> nx.MultiGraph:
        assignments = self._load_assignments(binning_tsv)
        if not assignments:
            logger.warning("Binning miner: no assignments loaded from %s", binning_tsv)
            return self.graph

        self._ensure_node_defaults()
        seeded_nodes = self._apply_input_assignments(assignments)
        rescued_nodes = self._rescue_unbinned_contigs(seeded_nodes)
        membership_edges = self._inject_bin_membership_edges()

        logger.info(
            "Binning miner complete: %d seeded nodes, %d rescued nodes, %d bin_membership edges.",
            len(seeded_nodes),
            rescued_nodes,
            membership_edges,
        )
        return self.graph

    def _load_assignments(self, binning_tsv: str | Path) -> dict[str, str]:
        path = Path(binning_tsv)
        if not path.exists():
            raise FileNotFoundError(f"Binning table not found: {path}")

        assignments: dict[str, str] = {}
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"Binning table has no header: {path}")

            if "Contig_ID" not in reader.fieldnames or "Bin_ID" not in reader.fieldnames:
                raise ValueError(
                    "Binning table must contain 'Contig_ID' and 'Bin_ID' columns. "
                    f"Found: {reader.fieldnames}"
                )

            for row in reader:
                contig_id = str(row.get("Contig_ID", "")).strip()
                bin_id = str(row.get("Bin_ID", "")).strip()
                if not contig_id or not bin_id:
                    continue
                assignments[contig_id] = bin_id

        return assignments

    def _ensure_node_defaults(self) -> None:
        for _, attrs in self.graph.nodes(data=True):
            attrs.setdefault("bin_id", None)
            attrs.setdefault("bin_status", "unbinned")
            attrs.setdefault("bin_source", "none")
            attrs.setdefault("rescued", False)
            attrs.setdefault("rescue_via", [])
            attrs.setdefault("rescue_support_count", 0)

    def _apply_input_assignments(self, assignments: dict[str, str]) -> set[str]:
        seeded_nodes: set[str] = set()
        for contig_id, bin_id in assignments.items():
            if not self.graph.has_node(contig_id):
                self.graph.add_node(contig_id, node_type="unknown", length=0)
            attrs = self.graph.nodes[contig_id]
            attrs["bin_id"] = bin_id
            attrs["bin_status"] = "binned"
            attrs["bin_source"] = "input"
            attrs["rescued"] = False
            attrs["rescue_via"] = []
            attrs["rescue_support_count"] = 0
            seeded_nodes.add(contig_id)
        return seeded_nodes

    def _rescue_unbinned_contigs(self, seeded_nodes: set[str]) -> int:
        rescued = 0
        for node_id in list(self.graph.nodes()):
            attrs = self.graph.nodes[node_id]
            if attrs.get("bin_id") not in (None, "", "NA"):
                continue

            support_by_bin: dict[str, set[str]] = defaultdict(set)
            support_types_by_bin: dict[str, set[str]] = defaultdict(set)
            for _, neighbor, edge_data in self.graph.edges(node_id, data=True):
                edge_type = str(edge_data.get("type", ""))
                if edge_type not in self.strong_evidence_types:
                    continue
                if neighbor not in seeded_nodes:
                    continue
                neighbor_bin = self.graph.nodes[neighbor].get("bin_id")
                if not neighbor_bin:
                    continue
                support_by_bin[str(neighbor_bin)].add(str(neighbor))
                support_types_by_bin[str(neighbor_bin)].add(edge_type)

            if not support_by_bin:
                continue

            score_rows = sorted(
                (
                    (bin_id, len(neighbors), len(support_types_by_bin.get(bin_id, set())))
                    for bin_id, neighbors in support_by_bin.items()
                ),
                key=lambda row: (row[1], row[2], row[0]),
                reverse=True,
            )
            best_bin, support_count, _ = score_rows[0]
            if support_count < self.min_rescue_support:
                continue
            if len(score_rows) > 1 and score_rows[1][1] == support_count:
                attrs["bin_conflict"] = [row[0] for row in score_rows if row[1] == support_count]
                continue

            attrs["bin_id"] = best_bin
            attrs["bin_status"] = "rescued"
            attrs["bin_source"] = "rescued"
            attrs["rescued"] = True
            attrs["rescue_via"] = sorted(support_types_by_bin[best_bin])
            attrs["rescue_support_count"] = support_count
            rescued += 1

        return rescued

    def _inject_bin_membership_edges(self) -> int:
        members_by_bin: dict[str, list[str]] = defaultdict(list)
        for node_id, attrs in self.graph.nodes(data=True):
            bin_id = attrs.get("bin_id")
            if bin_id:
                members_by_bin[str(bin_id)].append(str(node_id))

        added = 0
        for bin_id, members in members_by_bin.items():
            if len(members) < 2:
                continue

            if self.dense_membership:
                for idx, source in enumerate(members):
                    for target in members[idx + 1 :]:
                        if self._has_bin_membership(source, target, bin_id):
                            continue
                        self.graph.add_edge(
                            source,
                            target,
                            type="bin_membership",
                            bin_id=bin_id,
                            assignment_source=self._assignment_source(source, target),
                            rescue_support_count=max(
                                int(self.graph.nodes[source].get("rescue_support_count", 0)),
                                int(self.graph.nodes[target].get("rescue_support_count", 0)),
                            ),
                            source="binning_miner",
                        )
                        added += 1
                continue

            representative = max(members, key=self._node_length)
            for member in members:
                if member == representative:
                    continue
                if self._has_bin_membership(representative, member, bin_id):
                    continue
                self.graph.add_edge(
                    representative,
                    member,
                    type="bin_membership",
                    bin_id=bin_id,
                    assignment_source=self._assignment_source(representative, member),
                    rescue_support_count=max(
                        int(self.graph.nodes[representative].get("rescue_support_count", 0)),
                        int(self.graph.nodes[member].get("rescue_support_count", 0)),
                    ),
                    source="binning_miner",
                )
                added += 1

        return added

    def _assignment_source(self, source: str, target: str) -> str:
        source_kind = self.graph.nodes[source].get("bin_source", "none")
        target_kind = self.graph.nodes[target].get("bin_source", "none")
        if source_kind == "rescued" or target_kind == "rescued":
            return "rescued"
        return "input"

    def _has_bin_membership(self, u: str, v: str, bin_id: str) -> bool:
        edge_dict = self.graph.get_edge_data(u, v, default={})
        for edge_data in edge_dict.values():
            if edge_data.get("type") == "bin_membership" and str(edge_data.get("bin_id")) == str(bin_id):
                return True
        return False

    def _node_length(self, node_id: str) -> int:
        return int(self.graph.nodes[node_id].get("length", 0) or 0)
