"""
Unit tests for Module 3: Graph Exporter (TSV + HTML)
"""
import textwrap
from pathlib import Path

import networkx as nx
import pytest

from contigweaver.modules.graph_exporter import GraphExporter, export_graph


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


def test_export_graph_convenience(sample_graph, tmp_path):
    tsv, html = export_graph(
        sample_graph,
        tsv_path=tmp_path / "e.tsv",
        html_path=tmp_path / "n.html",
        viral_contigs={"virus_1"},
    )
    assert tsv.exists()
    assert html.exists()
