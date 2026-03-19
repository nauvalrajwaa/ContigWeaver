"""
Unit tests for Module 3: Graph Exporter (TSV + HTML)
"""
import textwrap
from pathlib import Path

import networkx as nx
import pytest

from contigweaver.modules.graph_exporter import (
    GraphExporter,
    export_graph,
    export_index_report,
)


@pytest.fixture()
def sample_graph() -> nx.MultiGraph:
    g = nx.MultiGraph()
    g.add_node("bact_1", length=1000, node_type="unknown")
    g.add_node("bact_2", length=2500, node_type="unknown")
    g.add_node("virus_1", length=500, node_type="viral")

    g.add_edge("bact_1", "bact_2",
               type="physical_overlap",
               overlap_cigar="55M",
               from_orient="+",
               to_orient="+")
    g.add_edge("bact_1", "virus_1",
               type="crispr_targeting",
               identity=97.5,
               coverage=0.95,
               evalue=1e-10)
    g.add_edge("bact_1", "bact_2",
               type="co_abundance_guild",
               weight=0.92,
               p_value=0.002,
               metabolic_match=True)
    return g


# ---------------------------------------------------------------------------
# TSV export tests
# ---------------------------------------------------------------------------

def test_tsv_header(sample_graph, tmp_path):
    exporter = GraphExporter(sample_graph)
    out = tmp_path / "edges.tsv"
    exporter.export_tsv(out)
    lines = out.read_text().strip().split("\n")
    header = lines[0]
    assert "Source_Contig" in header
    assert "Target_Contig" in header
    assert "Evidence_Type" in header
    assert "Weight_or_Identity" in header
    assert "Source_Length" in header
    assert "Target_Length" in header


def test_tsv_row_count(sample_graph, tmp_path):
    exporter = GraphExporter(sample_graph)
    out = tmp_path / "edges.tsv"
    exporter.export_tsv(out)
    lines = out.read_text().strip().split("\n")
    # 1 header + 3 edges
    assert len(lines) == 4


def test_tsv_evidence_type_mapping(sample_graph, tmp_path):
    exporter = GraphExporter(sample_graph)
    out = tmp_path / "edges.tsv"
    exporter.export_tsv(out)
    content = out.read_text()
    assert "Physical" in content
    assert "CRISPR" in content
    assert "Co-Abundance" in content


def test_tsv_evidence_type_mapping_includes_new_layers(tmp_path):
    graph = nx.MultiGraph()
    graph.add_node("a", length=800, node_type="unknown")
    graph.add_node("b", length=950, node_type="unknown")
    graph.add_node("c", length=1100, node_type="unknown")
    graph.add_edge("a", "b", type="taxonomic_match", taxonomy_rank="genus", taxonomy_label="Bacillus", confidence_min=0.97)
    graph.add_edge("a", "c", type="functional_operon", cas_genes="Cas9", support_mode="nearby", distance_hops=1, score=0.5)
    graph.add_edge("b", "c", type="bin_membership", bin_id="bin_1", assignment_source="rescued", rescue_support_count=2)

    out = tmp_path / "edges.tsv"
    GraphExporter(graph).export_tsv(out)
    content = out.read_text()
    assert "Taxonomic-Match" in content
    assert "Functional-Operon" in content
    assert "Bin-Membership" in content


# ---------------------------------------------------------------------------
# HTML export tests (requires pyvis)
# ---------------------------------------------------------------------------

pytest.importorskip("pyvis", reason="pyvis not installed")


def test_html_export_creates_file(sample_graph, tmp_path):
    exporter = GraphExporter(sample_graph)
    out = tmp_path / "network.html"
    exporter.export_html(out, viral_contigs={"virus_1"})
    assert out.exists()
    assert out.stat().st_size > 1000  # should be a real HTML file


def test_html_export_contains_node_labels(sample_graph, tmp_path):
    exporter = GraphExporter(sample_graph)
    out = tmp_path / "network.html"
    exporter.export_html(out, viral_contigs={"virus_1"})
    html = out.read_text()
    assert "bact_1" in html
    assert "virus_1" in html


def test_html_export_includes_edge_filter_controls_and_node_metadata(tmp_path):
    graph = nx.MultiGraph()
    graph.add_node(
        "host_1",
        length=1200,
        node_type="unknown",
        bin_id="bin_alpha",
        bin_status="rescued",
        rescued=True,
        taxonomy_label="Bacillus",
        taxonomy_rank="genus",
        taxonomy_confidence=0.97,
        cas_genes=["Cas9"],
        has_crispr=True,
        spacer_count=2,
    )
    graph.add_node("host_2", length=1500, node_type="unknown", bin_id="bin_alpha")
    graph.add_node("virus_1", length=600, node_type="viral")
    graph.add_edge("host_1", "virus_1", type="crispr_targeting", identity=99.0, coverage=1.0)
    graph.add_edge("host_1", "host_2", type="functional_operon", cas_genes="Cas9", support_mode="nearby", distance_hops=1, score=0.5)
    graph.add_edge("host_1", "host_2", type="taxonomic_match", taxonomy_rank="genus", taxonomy_label="Bacillus", confidence_min=0.95)
    graph.add_edge("host_1", "host_2", type="bin_membership", bin_id="bin_alpha", assignment_source="rescued", rescue_support_count=2)

    out = tmp_path / "network.html"
    GraphExporter(graph).export_html(out, viral_contigs={"virus_1"})
    html = out.read_text()
    assert "Edge Filters" in html
    assert 'data-edge-type="taxonomic_match"' in html
    assert 'data-edge-type="functional_operon"' in html
    assert 'data-edge-type="bin_membership"' in html
    assert "Bin: bin_alpha" in html
    assert "Taxonomy: Bacillus" in html


def test_focus_nodes_ignores_bin_membership_only_edges():
    graph = nx.MultiGraph()
    graph.add_node("a", length=500, node_type="unknown")
    graph.add_node("b", length=520, node_type="unknown")
    graph.add_edge("a", "b", type="bin_membership", bin_id="bin_1")
    exporter = GraphExporter(graph)
    assert exporter._focus_nodes(graph) == []


def test_export_graph_convenience(sample_graph, tmp_path):
    tsv, html = export_graph(
        sample_graph,
        tsv_path=tmp_path / "e.tsv",
        html_path=tmp_path / "n.html",
        viral_contigs={"virus_1"},
    )
    assert tsv.exists()
    assert html.exists()


def test_export_index_report_creates_sections(sample_graph, tmp_path):
    network_html = tmp_path / "network.html"
    edges_tsv = tmp_path / "edges.tsv"
    GraphExporter(sample_graph).export_html(network_html, viral_contigs={"virus_1"})
    GraphExporter(sample_graph).export_tsv(edges_tsv)

    related_html = tmp_path / "quast.html"
    related_html.write_text("<html><body><h1>QUAST</h1></body></html>")

    index_html = export_index_report(
        sample_graph,
        output_path=tmp_path / "index.html",
        network_html_path=network_html,
        edges_tsv_path=edges_tsv,
        viral_contigs={"virus_1"},
        input_paths={"Assembly graph": tmp_path / "assembly.gfa"},
        related_html_paths=[related_html],
    )

    content = index_html.read_text()
    assert index_html.exists()
    assert "ContigWeaver Comprehensive Report" in content
    assert "Interactive Network" in content
    assert "Related HTML Reports" in content
    assert "quast.html" in content


def test_prepare_html_graph_samples_large_graph_around_focus_nodes():
    graph = nx.MultiGraph()
    for idx in range(60):
        node_id = str(idx)
        graph.add_node(node_id, node_type="unknown", length=500 + idx)
        if idx > 0:
            graph.add_edge(str(idx - 1), node_id, type="physical_overlap", overlap_cigar="55M")

    graph.add_node("NODE_host", node_type="unknown", length=1800)
    graph.add_node("NODE_virus", node_type="viral", length=2200)
    graph.add_edge("NODE_host", "NODE_virus", type="crispr_targeting", identity=98.0, coverage=1.0)
    graph.add_edge("NODE_host", "10", type="segment_membership", anchor_length=300, match_type="substring", orientation="+")
    graph.add_edge("NODE_virus", "11", type="segment_membership", anchor_length=280, match_type="substring", orientation="+")

    exporter = GraphExporter(graph)
    html_graph, meta = exporter._prepare_html_graph(max_nodes=12, max_edges=16, focus_hops=1)

    assert meta["sampled"] is False
    assert meta["hairball_filtered"] is True
    assert meta["ego_center"] == "NODE_virus"
    assert html_graph.number_of_nodes() <= 12
    assert html_graph.number_of_edges() <= 16
    assert html_graph.has_node("NODE_host")
    assert html_graph.has_node("NODE_virus")
    edge_types = {data["type"] for _, _, data in html_graph.edges(data=True)}
    assert "crispr_targeting" in edge_types
    assert "segment_membership" in edge_types


def test_filter_physical_hairball_keeps_only_host_virus_physical_edges():
    graph = nx.MultiGraph()
    graph.add_node("host_a", node_type="unknown", length=1000)
    graph.add_node("host_b", node_type="unknown", length=1100)
    graph.add_node("virus_a", node_type="viral", length=900)
    graph.add_node("virus_b", node_type="viral", length=950)
    graph.add_edge("host_a", "host_b", type="physical_overlap", overlap_cigar="55M")
    graph.add_edge("virus_a", "virus_b", type="physical_overlap", overlap_cigar="40M")
    graph.add_edge("host_a", "virus_a", type="physical_overlap", overlap_cigar="25M")

    filtered = GraphExporter(graph)._filter_physical_hairball({"virus_a", "virus_b"})
    physical_pairs = {
        frozenset((str(u), str(v)))
        for u, v, data in filtered.edges(data=True)
        if data.get("type") == "physical_overlap"
    }

    assert frozenset(("host_a", "virus_a")) in physical_pairs
    assert frozenset(("host_a", "host_b")) not in physical_pairs
    assert frozenset(("virus_a", "virus_b")) not in physical_pairs


def test_html_export_filter_defaults_and_auto_summary(tmp_path):
    graph = nx.MultiGraph()
    graph.add_node("host_1", node_type="unknown", length=1200)
    graph.add_node("host_2", node_type="unknown", length=1400)
    graph.add_node("virus_1", node_type="viral", length=800)
    graph.add_edge("host_1", "virus_1", type="crispr_targeting", identity=98.5, coverage=0.95)
    graph.add_edge("host_1", "host_2", type="co_abundance_guild", weight=0.9)
    graph.add_edge("host_1", "host_2", type="physical_overlap", overlap_cigar="80M")
    graph.add_edge("host_1", "virus_1", type="physical_overlap", overlap_cigar="45M")
    graph.add_edge("virus_1", "host_2", type="segment_membership", anchor_length=210)

    out = tmp_path / "network.html"
    GraphExporter(graph).export_html(out, viral_contigs={"virus_1"})
    html = out.read_text()

    assert 'data-edge-type="crispr_targeting" checked' in html
    assert 'data-edge-type="segment_membership" checked' in html
    assert 'data-edge-type="co_abundance_guild" checked' in html
    assert 'data-edge-type="physical_overlap">' in html
    assert 'data-edge-type="physical_overlap" checked' not in html
    assert "Ringkasan Ekosistem" in html


def test_export_index_report_includes_multi_evidence_key_findings(tmp_path):
    graph = nx.MultiGraph()
    graph.add_node("host_1", node_type="unknown", length=1500)
    graph.add_node("host_2", node_type="unknown", length=1300)
    graph.add_node("virus_1", node_type="viral", length=900)
    graph.add_edge("host_1", "virus_1", type="crispr_targeting", identity=99.0, coverage=0.98)
    graph.add_edge("host_1", "virus_1", type="co_abundance_guild", weight=0.91)
    graph.add_edge("host_2", "virus_1", type="segment_membership", anchor_length=250)

    network_html = tmp_path / "network.html"
    edges_tsv = tmp_path / "edges.tsv"
    network_html.write_text("<html><body>network</body></html>")
    edges_tsv.write_text("Source_Contig\tTarget_Contig\tEvidence_Type\n")

    index_html = export_index_report(
        graph,
        output_path=tmp_path / "index.html",
        network_html_path=network_html,
        edges_tsv_path=edges_tsv,
        viral_contigs={"virus_1"},
    )

    content = index_html.read_text()
    assert "Key Findings" in content
    assert "host_1" in content
    assert "virus_1" in content
    assert "crispr_targeting" in content
    assert "co_abundance_guild" in content
