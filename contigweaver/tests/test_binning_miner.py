from pathlib import Path

import networkx as nx
import pytest

from contigweaver.modules.binning_miner import BinningMiner


def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> Path:
    body = ["\t".join(header)] + ["\t".join(row) for row in rows]
    path.write_text("\n".join(body) + "\n")
    return path


def _edge_count(graph: nx.MultiGraph, edge_type: str) -> int:
    return sum(1 for _, _, data in graph.edges(data=True) if data.get("type") == edge_type)


def test_binning_miner_rescues_unbinned_contigs_and_adds_membership(tmp_path: Path):
    graph = nx.MultiGraph()
    graph.add_node("contig_a", length=900)
    graph.add_node("contig_b", length=1200)
    graph.add_node("contig_c", length=1800)

    graph.add_edge("contig_a", "contig_b", type="physical_overlap", overlap_cigar="50M")
    graph.add_edge("contig_a", "contig_c", type="co_abundance_guild", weight=0.91, p_value=0.01)

    bins = _write_tsv(
        tmp_path / "bins.tsv",
        ["Contig_ID", "Bin_ID"],
        [["contig_a", "bin_001"], ["contig_b", "bin_001"]],
    )

    miner = BinningMiner(graph=graph)
    miner.run(bins)

    attrs_c = graph.nodes["contig_c"]
    assert attrs_c["bin_id"] == "bin_001"
    assert attrs_c["bin_status"] == "rescued"
    assert attrs_c["rescued"] is True
    assert "co_abundance_guild" in attrs_c["rescue_via"]

    assert _edge_count(graph, "bin_membership") >= 2


def test_binning_miner_leaves_conflicted_node_unbinned(tmp_path: Path):
    graph = nx.MultiGraph()
    graph.add_node("contig_a", length=900)
    graph.add_node("contig_b", length=950)
    graph.add_node("contig_d", length=700)

    graph.add_edge("contig_d", "contig_a", type="physical_overlap", overlap_cigar="60M")
    graph.add_edge("contig_d", "contig_b", type="physical_overlap", overlap_cigar="61M")

    bins = _write_tsv(
        tmp_path / "bins.tsv",
        ["Contig_ID", "Bin_ID"],
        [["contig_a", "bin_A"], ["contig_b", "bin_B"]],
    )

    miner = BinningMiner(graph=graph)
    miner.run(bins)

    attrs_d = graph.nodes["contig_d"]
    assert attrs_d.get("bin_id") in (None, "")
    assert set(attrs_d.get("bin_conflict", [])) == {"bin_A", "bin_B"}


def test_binning_miner_requires_columns(tmp_path: Path):
    graph = nx.MultiGraph()
    bad_bins = _write_tsv(
        tmp_path / "bad_bins.tsv",
        ["Contig", "Bin"],
        [["contig_a", "bin_001"]],
    )

    miner = BinningMiner(graph=graph)
    with pytest.raises(ValueError, match="Contig_ID"):
        miner.run(bad_bins)
