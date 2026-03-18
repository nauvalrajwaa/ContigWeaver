from pathlib import Path

import networkx as nx
import pytest

from contigweaver.modules.annotation_miner import AnnotationMiner


def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> Path:
    body = ["\t".join(header)] + ["\t".join(row) for row in rows]
    path.write_text("\n".join(body) + "\n")
    return path


def _edge_count(graph: nx.MultiGraph, edge_type: str) -> int:
    return sum(1 for _, _, data in graph.edges(data=True) if data.get("type") == edge_type)


def test_annotation_miner_adds_metadata_and_edges(tmp_path: Path):
    graph = nx.MultiGraph()
    graph.add_node("contig_a", length=850, node_type="unknown")
    graph.add_node("contig_b", length=1400, node_type="unknown")
    graph.add_node("contig_c", length=1800, node_type="unknown")
    graph.add_node("contig_d", length=2100, node_type="unknown")
    graph.add_edge("contig_a", "contig_b", type="physical_overlap", overlap_cigar="55M")

    annotations = _write_tsv(
        tmp_path / "annotations.tsv",
        [
            "Contig_ID",
            "taxonomy_label",
            "taxonomy_rank",
            "taxonomy_confidence",
            "functional_terms",
            "cas_genes",
            "has_crispr",
            "spacer_count",
        ],
        [
            ["contig_a", "Bacillus", "genus", "0.99", "CRISPR repeat", "", "true", "2"],
            ["contig_b", "Bacillus", "genus", "0.98", "DNA repair", "Cas9", "false", "0"],
            ["contig_c", "Bacillus", "genus", "0.97", "Endonuclease", "Cas12", "false", "0"],
            ["contig_d", "Bacillus", "genus", "0.40", "Other", "", "false", "0"],
        ],
    )

    miner = AnnotationMiner(graph=graph, taxonomy_confidence_threshold=0.90, nearby_hops=2)
    miner.run(annotations)

    attrs_a = graph.nodes["contig_a"]
    assert attrs_a["taxonomy_label"] == "Bacillus"
    assert attrs_a["has_crispr"] is True
    assert attrs_a["spacer_count"] == 2

    attrs_b = graph.nodes["contig_b"]
    assert attrs_b["cas_genes"] == ["Cas9"]

    assert _edge_count(graph, "functional_operon") >= 2
    assert _edge_count(graph, "taxonomic_match") == 1


def test_annotation_miner_requires_contig_id_column(tmp_path: Path):
    graph = nx.MultiGraph()
    path = _write_tsv(
        tmp_path / "bad_annotations.tsv",
        ["Wrong_ID", "taxonomy_label"],
        [["x", "Bacillus"]],
    )

    miner = AnnotationMiner(graph=graph)
    with pytest.raises(ValueError, match="Contig_ID"):
        miner.run(path)
