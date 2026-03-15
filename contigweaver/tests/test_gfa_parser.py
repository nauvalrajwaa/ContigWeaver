"""
Unit tests for Module 1: GFA Parser
"""
import tempfile
from pathlib import Path

import networkx as nx
import pytest

from contigweaver.modules.gfa_parser import GFAParser, parse_gfa


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

GFA_CONTENT = """\
H\tVN:Z:1.0
S\tcontig_1\t*\tLN:i:500
S\tcontig_2\t*\tLN:i:1200
S\tcontig_3\t*\tLN:i:800
L\tcontig_1\t+\tcontig_2\t+\t55M
L\tcontig_2\t-\tcontig_3\t+\t30M
"""


@pytest.fixture()
def gfa_file(tmp_path: Path) -> Path:
    p = tmp_path / "test.gfa"
    p.write_text(GFA_CONTENT)
    return p


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_parse_segments(gfa_file: Path):
    parser = GFAParser()
    graph = parser.parse(gfa_file)

    assert graph.has_node("contig_1")
    assert graph.has_node("contig_2")
    assert graph.has_node("contig_3")
    assert graph.nodes["contig_1"]["length"] == 500
    assert graph.nodes["contig_2"]["length"] == 1200
    assert graph.nodes["contig_3"]["length"] == 800


def test_parse_links(gfa_file: Path):
    parser = GFAParser()
    graph = parser.parse(gfa_file)

    assert graph.number_of_edges() == 2
    edge_types = {
        d["type"]
        for _, _, d in graph.edges(data=True)
    }
    assert edge_types == {"physical_overlap"}


def test_parse_link_attributes(gfa_file: Path):
    parser = GFAParser()
    graph = parser.parse(gfa_file)

    # Find the edge between contig_1 and contig_2
    edges = graph.get_edge_data("contig_1", "contig_2")
    assert edges is not None
    # MultiGraph: get_edge_data returns dict of {key: attrs}
    first_edge = list(edges.values())[0]
    assert first_edge["overlap_cigar"] == "55M"
    assert first_edge["from_orient"] == "+"
    assert first_edge["to_orient"] == "+"


def test_file_not_found_raises():
    parser = GFAParser()
    with pytest.raises(FileNotFoundError):
        parser.parse("/non/existent/file.gfa")


def test_empty_gfa_raises_value_error(tmp_path: Path):
    empty = tmp_path / "empty.gfa"
    empty.write_text("H\tVN:Z:1.0\n")  # header only, no S lines
    parser = GFAParser()
    with pytest.raises(ValueError, match="No Segment"):
        parser.parse(empty)


def test_sequence_based_length(tmp_path: Path):
    """When there is no LN tag, length should be inferred from the sequence."""
    gfa = tmp_path / "seq.gfa"
    gfa.write_text("S\tnode_a\tACGTACGT\n")
    parser = GFAParser()
    graph = parser.parse(gfa)
    assert graph.nodes["node_a"]["length"] == 8


def test_inject_into_existing_graph(gfa_file: Path):
    """Parsing should extend an existing graph."""
    existing = nx.MultiGraph()
    existing.add_node("pre_existing", length=999, node_type="bacterial")
    parser = GFAParser(graph=existing)
    graph = parser.parse(gfa_file)
    assert graph.has_node("pre_existing")
    assert graph.has_node("contig_1")


def test_convenience_function(gfa_file: Path):
    graph = parse_gfa(gfa_file)
    assert isinstance(graph, nx.MultiGraph)
    assert graph.number_of_nodes() == 3
