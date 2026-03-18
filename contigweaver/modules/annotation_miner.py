from __future__ import annotations

import csv
import logging
import re
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

import networkx as nx

logger = logging.getLogger(__name__)

TAXONOMY_LABEL_COLUMNS = (
    "taxonomy_label",
    "taxonomy",
    "taxon",
    "classification",
    "kraken_taxonomy",
)
TAXONOMY_RANK_COLUMNS = ("taxonomy_rank", "rank", "tax_rank")
TAXONOMY_CONFIDENCE_COLUMNS = (
    "taxonomy_confidence",
    "tax_confidence",
    "confidence",
    "score",
)
FUNCTIONAL_TERMS_COLUMNS = ("functional_terms", "function", "gene_terms", "genes")
CAS_GENE_COLUMNS = ("cas_genes", "cas_terms")
CRISPR_FLAG_COLUMNS = ("has_crispr", "crispr_positive")
SPACER_COUNT_COLUMNS = ("spacer_count", "crispr_spacer_count", "spacers")

CAS_GENE_PATTERN = re.compile(r"\bCAS(?:\d+|[A-Z])?\b", re.IGNORECASE)


@dataclass(slots=True)
class AnnotationRecord:
    contig_id: str
    taxonomy_label: str | None
    taxonomy_rank: str | None
    taxonomy_confidence: float
    functional_terms: list[str]
    cas_genes: list[str]
    has_crispr: bool
    spacer_count: int


class AnnotationMiner:
    def __init__(
        self,
        graph: nx.MultiGraph,
        taxonomy_confidence_threshold: float = 0.90,
        small_contig_max_bp: int = 5000,
        nearby_hops: int = 2,
        max_taxonomic_edges: int = 5000,
        max_operon_edges: int = 5000,
    ):
        self.graph = graph
        self.taxonomy_confidence_threshold = taxonomy_confidence_threshold
        self.small_contig_max_bp = small_contig_max_bp
        self.nearby_hops = nearby_hops
        self.max_taxonomic_edges = max_taxonomic_edges
        self.max_operon_edges = max_operon_edges

    def run(self, annotations_tsv: str | Path) -> nx.MultiGraph:
        records = self._load_records(annotations_tsv)
        if not records:
            logger.warning("Annotation miner: no records loaded from %s", annotations_tsv)
            return self.graph

        self._apply_node_metadata(records)
        taxonomic_added = self._inject_taxonomic_matches(records)
        functional_added = self._inject_functional_operons(records)

        logger.info(
            "Annotation miner complete: %d records, %d taxonomic_match edges, %d functional_operon edges.",
            len(records),
            taxonomic_added,
            functional_added,
        )
        return self.graph

    def _load_records(self, annotations_tsv: str | Path) -> dict[str, AnnotationRecord]:
        path = Path(annotations_tsv)
        if not path.exists():
            raise FileNotFoundError(f"Annotation data not found: {path}")

        records: dict[str, AnnotationRecord] = {}
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"Annotation table has no header: {path}")
            if "Contig_ID" not in reader.fieldnames:
                raise ValueError(
                    "Annotation table must contain a 'Contig_ID' column. "
                    f"Found: {reader.fieldnames}"
                )

            for row in reader:
                contig_id = str(row.get("Contig_ID", "")).strip()
                if not contig_id:
                    continue

                functional_terms = self._split_terms(
                    self._coalesce(row, FUNCTIONAL_TERMS_COLUMNS)
                )
                explicit_cas = self._split_terms(self._coalesce(row, CAS_GENE_COLUMNS))
                detected_cas = sorted(
                    {
                        term
                        for term in functional_terms
                        if CAS_GENE_PATTERN.search(term) is not None
                    }
                )
                cas_genes = sorted(set(explicit_cas) | set(detected_cas))

                has_crispr = self._parse_bool(self._coalesce(row, CRISPR_FLAG_COLUMNS))
                spacer_count = self._parse_int(self._coalesce(row, SPACER_COUNT_COLUMNS), default=0)

                record = AnnotationRecord(
                    contig_id=contig_id,
                    taxonomy_label=self._clean_text(self._coalesce(row, TAXONOMY_LABEL_COLUMNS)),
                    taxonomy_rank=self._clean_text(self._coalesce(row, TAXONOMY_RANK_COLUMNS)),
                    taxonomy_confidence=self._parse_float(
                        self._coalesce(row, TAXONOMY_CONFIDENCE_COLUMNS),
                        default=0.0,
                    ),
                    functional_terms=functional_terms,
                    cas_genes=cas_genes,
                    has_crispr=has_crispr,
                    spacer_count=spacer_count,
                )
                records[contig_id] = record

        return records

    def _apply_node_metadata(self, records: dict[str, AnnotationRecord]) -> None:
        for contig_id, record in records.items():
            if not self.graph.has_node(contig_id):
                self.graph.add_node(contig_id, node_type="unknown", length=0)

            attrs = self.graph.nodes[contig_id]
            if record.taxonomy_label is not None:
                attrs["taxonomy_label"] = record.taxonomy_label
            if record.taxonomy_rank is not None:
                attrs["taxonomy_rank"] = record.taxonomy_rank
            attrs["taxonomy_confidence"] = record.taxonomy_confidence
            attrs["functional_terms"] = record.functional_terms
            attrs["cas_genes"] = record.cas_genes

            inferred_crispr = record.has_crispr or record.spacer_count > 0 or self._has_crispr_edge(contig_id)
            attrs["has_crispr"] = inferred_crispr
            attrs["spacer_count"] = max(int(attrs.get("spacer_count", 0)), record.spacer_count)

    def _inject_taxonomic_matches(self, records: dict[str, AnnotationRecord]) -> int:
        component_id: dict[str, int] = {}
        for idx, component in enumerate(nx.connected_components(nx.Graph(self.graph))):
            for node in component:
                component_id[str(node)] = idx

        groups: dict[tuple[str, str], list[AnnotationRecord]] = defaultdict(list)
        for record in records.values():
            if (
                record.taxonomy_label is None
                or record.taxonomy_rank is None
                or record.taxonomy_confidence < self.taxonomy_confidence_threshold
            ):
                continue
            if record.contig_id not in component_id:
                continue
            groups[(record.taxonomy_rank, record.taxonomy_label)].append(record)

        edges_added = 0
        for (rank, label), members in groups.items():
            by_component: dict[int, list[AnnotationRecord]] = defaultdict(list)
            for record in members:
                cid = component_id.get(record.contig_id)
                if cid is None:
                    continue
                by_component[cid].append(record)
            if len(by_component) < 2:
                continue

            representative_by_component: dict[int, AnnotationRecord] = {}
            for cid, cid_records in by_component.items():
                representative_by_component[cid] = max(
                    cid_records,
                    key=lambda rec: (rec.taxonomy_confidence, self._node_length(rec.contig_id)),
                )

            for cid_a, cid_b in combinations(sorted(representative_by_component), 2):
                if edges_added >= self.max_taxonomic_edges:
                    logger.warning(
                        "Taxonomic edge cap (%d) reached; truncating additional matches.",
                        self.max_taxonomic_edges,
                    )
                    return edges_added

                rec_a = representative_by_component[cid_a]
                rec_b = representative_by_component[cid_b]
                if self._has_edge_between(rec_a.contig_id, rec_b.contig_id, "taxonomic_match"):
                    continue

                self.graph.add_edge(
                    rec_a.contig_id,
                    rec_b.contig_id,
                    type="taxonomic_match",
                    taxonomy_rank=rank,
                    taxonomy_label=label,
                    confidence_a=rec_a.taxonomy_confidence,
                    confidence_b=rec_b.taxonomy_confidence,
                    confidence_min=min(rec_a.taxonomy_confidence, rec_b.taxonomy_confidence),
                    source="annotation_miner",
                )
                edges_added += 1

        return edges_added

    def _inject_functional_operons(self, records: dict[str, AnnotationRecord]) -> int:
        simple_graph = nx.Graph(self.graph)
        component_id: dict[str, int] = {}
        for idx, component in enumerate(nx.connected_components(simple_graph)):
            for node in component:
                component_id[str(node)] = idx

        crispr_nodes = [
            record
            for record in records.values()
            if self._is_small_crispr_contig(record)
        ]
        cas_nodes = [
            record
            for record in records.values()
            if len(record.cas_genes) > 0
        ]

        edges_added = 0
        for crispr_record in crispr_nodes:
            for cas_record in cas_nodes:
                if crispr_record.contig_id == cas_record.contig_id:
                    continue
                if edges_added >= self.max_operon_edges:
                    logger.warning(
                        "Functional operon edge cap (%d) reached; truncating additional links.",
                        self.max_operon_edges,
                    )
                    return edges_added
                if self._has_edge_between(
                    crispr_record.contig_id,
                    cas_record.contig_id,
                    "functional_operon",
                ):
                    continue

                support_mode, distance_hops = self._operon_support_mode(
                    crispr_record,
                    cas_record,
                    simple_graph,
                    component_id,
                )
                if support_mode is None:
                    continue

                score = 1.0 / (1.0 + float(distance_hops)) if distance_hops >= 0 else 0.5
                self.graph.add_edge(
                    crispr_record.contig_id,
                    cas_record.contig_id,
                    type="functional_operon",
                    spacer_count=max(crispr_record.spacer_count, 1),
                    cas_genes=",".join(cas_record.cas_genes),
                    support_mode=support_mode,
                    distance_hops=distance_hops,
                    score=score,
                    source="annotation_miner",
                )
                edges_added += 1

        return edges_added

    def _operon_support_mode(
        self,
        crispr_record: AnnotationRecord,
        cas_record: AnnotationRecord,
        simple_graph: nx.Graph,
        component_id: dict[str, int],
    ) -> tuple[str | None, int]:
        a = crispr_record.contig_id
        b = cas_record.contig_id
        cid_a = component_id.get(a)
        cid_b = component_id.get(b)
        if cid_a is None or cid_b is None:
            return None, -1

        if cid_a == cid_b:
            try:
                distance = nx.shortest_path_length(simple_graph, source=a, target=b)
            except nx.NetworkXNoPath:
                return None, -1
            if distance <= self.nearby_hops:
                return "nearby", int(distance)
            return None, -1

        if self._taxonomy_match(crispr_record, cas_record):
            return "disconnected", -1
        return None, -1

    def _taxonomy_match(self, a: AnnotationRecord, b: AnnotationRecord) -> bool:
        if (
            a.taxonomy_label is None
            or b.taxonomy_label is None
            or a.taxonomy_rank is None
            or b.taxonomy_rank is None
        ):
            return False
        if a.taxonomy_label != b.taxonomy_label or a.taxonomy_rank != b.taxonomy_rank:
            return False
        return (
            a.taxonomy_confidence >= self.taxonomy_confidence_threshold
            and b.taxonomy_confidence >= self.taxonomy_confidence_threshold
        )

    def _is_small_crispr_contig(self, record: AnnotationRecord) -> bool:
        has_crispr = record.has_crispr or record.spacer_count > 0 or self._has_crispr_edge(record.contig_id)
        return has_crispr and self._node_length(record.contig_id) <= self.small_contig_max_bp

    def _has_crispr_edge(self, node_id: str) -> bool:
        if not self.graph.has_node(node_id):
            return False
        for _, _, edge_data in self.graph.edges(node_id, data=True):
            if edge_data.get("type") == "crispr_targeting":
                return True
        return False

    def _node_length(self, node_id: str) -> int:
        if not self.graph.has_node(node_id):
            return 0
        return int(self.graph.nodes[node_id].get("length", 0) or 0)

    def _has_edge_between(self, u: str, v: str, edge_type: str) -> bool:
        edge_dict = self.graph.get_edge_data(u, v, default={})
        return any(data.get("type") == edge_type for data in edge_dict.values())

    @staticmethod
    def _coalesce(row: dict[str, str], candidates: Iterable[str]) -> str | None:
        for key in candidates:
            value = row.get(key)
            if value is not None and value.strip():
                return value.strip()
        return None

    @staticmethod
    def _split_terms(value: str | None) -> list[str]:
        if value is None:
            return []
        chunks = re.split(r"[,;|]", value)
        terms = [item.strip() for item in chunks if item.strip()]
        unique: list[str] = []
        seen: set[str] = set()
        for term in terms:
            term_upper = term.upper()
            if term_upper in seen:
                continue
            seen.add(term_upper)
            unique.append(term)
        return unique

    @staticmethod
    def _parse_bool(value: str | None) -> bool:
        if value is None:
            return False
        return value.strip().lower() in {"1", "true", "yes", "y", "t"}

    @staticmethod
    def _parse_int(value: str | None, default: int) -> int:
        if value is None:
            return default
        try:
            return int(float(value))
        except ValueError:
            return default

    @staticmethod
    def _parse_float(value: str | None, default: float) -> float:
        if value is None:
            return default
        try:
            return float(value)
        except ValueError:
            return default

    @staticmethod
    def _clean_text(value: str | None) -> str | None:
        if value is None:
            return None
        cleaned = value.strip()
        return cleaned if cleaned else None
