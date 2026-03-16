import gzip
import logging
from collections import defaultdict
from pathlib import Path
from typing import Sequence

import networkx as nx

logger = logging.getLogger(__name__)


def _open_text(path: str | Path):
    path = Path(path)
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1]


class ContigGraphReconciler:
    def __init__(
        self,
        graph: nx.MultiGraph,
        kmer_size: int = 31,
        max_anchors_per_contig: int = 8,
    ):
        self.graph = graph
        self.kmer_size = kmer_size
        self.max_anchors_per_contig = max_anchors_per_contig

    def reconcile(
        self,
        gfa_path: str | Path,
        fasta_paths: Sequence[str | Path],
    ) -> dict[str, int]:
        candidate_ids = self._candidate_contig_ids()
        if not candidate_ids:
            return {
                "candidate_contigs": 0,
                "resolved_contigs": 0,
                "segment_membership_edges": 0,
            }

        contig_sequences = self._load_contig_sequences(candidate_ids, fasta_paths)
        if not contig_sequences:
            logger.info("No candidate contig sequences found for reconciliation.")
            return {
                "candidate_contigs": len(candidate_ids),
                "resolved_contigs": 0,
                "segment_membership_edges": 0,
            }

        anchors = self._find_segment_anchors(gfa_path, contig_sequences)
        added_edges = 0
        resolved_contigs = 0

        for contig_id, matches in anchors.items():
            added_for_contig = 0
            for match in matches[: self.max_anchors_per_contig]:
                segment_id = str(match["segment_id"])
                if not self.graph.has_node(segment_id):
                    continue
                if self._has_membership_edge(contig_id, segment_id):
                    continue

                self.graph.add_edge(
                    contig_id,
                    segment_id,
                    type="segment_membership",
                    anchor_length=match["anchor_length"],
                    match_type=match["match_type"],
                    orientation=match["orientation"],
                )
                added_edges += 1
                added_for_contig += 1

            if added_for_contig:
                resolved_contigs += 1

        logger.info(
            "Contig reconciliation anchored %d/%d contigs with %d segment_membership edges.",
            resolved_contigs,
            len(candidate_ids),
            added_edges,
        )
        return {
            "candidate_contigs": len(candidate_ids),
            "resolved_contigs": resolved_contigs,
            "segment_membership_edges": added_edges,
        }

    def _candidate_contig_ids(self) -> list[str]:
        candidate_ids: set[str] = set()
        for node, attrs in self.graph.nodes(data=True):
            node_id = str(node)
            if not node_id.startswith("NODE_"):
                continue

            edge_types = {
                data.get("type")
                for _, _, data in self.graph.edges(node, data=True)
            }
            if attrs.get("node_type") == "viral" or any(
                edge_type and edge_type != "segment_membership"
                for edge_type in edge_types
            ):
                candidate_ids.add(node_id)

        return sorted(candidate_ids)

    def _load_contig_sequences(
        self,
        candidate_ids: list[str],
        fasta_paths: Sequence[str | Path],
    ) -> dict[str, str]:
        remaining = set(candidate_ids)
        sequences: dict[str, str] = {}

        for fasta_path in fasta_paths:
            fasta_path = Path(fasta_path)
            if not fasta_path.exists() or not remaining:
                continue

            current_id = None
            chunks: list[str] = []
            with _open_text(fasta_path) as fh:
                for raw_line in fh:
                    line = raw_line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if current_id in remaining:
                            sequences[current_id] = "".join(chunks).upper()
                            remaining.remove(current_id)
                        current_id = line[1:].split()[0]
                        chunks = []
                    else:
                        chunks.append(line)

                if current_id in remaining:
                    sequences[current_id] = "".join(chunks).upper()
                    remaining.remove(current_id)

        return sequences

    def _find_segment_anchors(
        self,
        gfa_path: str | Path,
        contig_sequences: dict[str, str],
    ) -> dict[str, list[dict[str, str | int]]]:
        if not contig_sequences:
            return {}

        kmer_index: dict[str, set[str]] = defaultdict(set)
        short_sequences: dict[str, set[str]] = defaultdict(set)

        for contig_id, sequence in contig_sequences.items():
            variants = {sequence, _reverse_complement(sequence)}
            for variant in variants:
                if len(variant) >= self.kmer_size:
                    for idx in range(len(variant) - self.kmer_size + 1):
                        kmer_index[variant[idx : idx + self.kmer_size]].add(contig_id)
                else:
                    short_sequences[variant].add(contig_id)

        anchors: dict[str, list[dict[str, str | int]]] = defaultdict(list)

        with _open_text(gfa_path) as fh:
            for raw_line in fh:
                if not raw_line.startswith("S\t"):
                    continue
                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue

                segment_id = fields[1]
                segment_seq = fields[2].upper()
                if segment_seq == "*":
                    continue

                candidate_ids: set[str] = set()
                variants = (
                    (segment_seq, "+"),
                    (_reverse_complement(segment_seq), "-"),
                )
                for variant_seq, _ in variants:
                    if len(variant_seq) >= self.kmer_size:
                        candidate_ids.update(kmer_index.get(variant_seq[: self.kmer_size], set()))
                    else:
                        candidate_ids.update(short_sequences.get(variant_seq, set()))

                if not candidate_ids:
                    continue

                for contig_id in candidate_ids:
                    contig_seq = contig_sequences[contig_id]
                    match = self._match_segment_to_contig(segment_seq, contig_seq)
                    if match is None:
                        continue

                    anchors[contig_id].append(
                        {
                            "segment_id": segment_id,
                            "anchor_length": len(segment_seq),
                            "match_type": match["match_type"],
                            "orientation": match["orientation"],
                        }
                    )

        resolved: dict[str, list[dict[str, str | int]]] = {}
        for contig_id, matches in anchors.items():
            by_segment: dict[str, dict[str, str | int]] = {}
            for match in sorted(matches, key=lambda item: int(item["anchor_length"]), reverse=True):
                segment_id = str(match["segment_id"])
                by_segment.setdefault(segment_id, match)
            resolved[contig_id] = list(by_segment.values())

        return resolved

    @staticmethod
    def _match_segment_to_contig(segment_seq: str, contig_seq: str) -> dict[str, str] | None:
        reverse_seq = _reverse_complement(segment_seq)
        if segment_seq == contig_seq:
            return {"match_type": "exact", "orientation": "+"}
        if reverse_seq == contig_seq:
            return {"match_type": "exact", "orientation": "-"}
        if segment_seq in contig_seq:
            return {"match_type": "substring", "orientation": "+"}
        if reverse_seq in contig_seq:
            return {"match_type": "substring", "orientation": "-"}
        return None

    def _has_membership_edge(self, contig_id: str, segment_id: str) -> bool:
        if not self.graph.has_edge(contig_id, segment_id):
            return False
        for edge_data in self.graph.get_edge_data(contig_id, segment_id).values():
            if edge_data.get("type") == "segment_membership":
                return True
        return False
