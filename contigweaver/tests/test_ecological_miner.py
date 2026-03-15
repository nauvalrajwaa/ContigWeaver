"""
Unit tests for Module 4/5/6: EcologicalMiner
"""
import textwrap
from pathlib import Path

import networkx as nx
import pytest

from contigweaver.modules.ecological_miner import (
    CoAbundanceCorrelator,
    EcologicalMiner,
    EcologicalGraphInjector,
    PathwayComplementarityChecker,
    SPEARMAN_THRESHOLD,
    P_VALUE_THRESHOLD,
)


# ---------------------------------------------------------------------------
# Coverage table fixtures
# ---------------------------------------------------------------------------

# 5 samples, 4 contigs — contig_1 & contig_2 are perfectly correlated
COVERAGE_TSV = textwrap.dedent("""\
    Contig_ID\tsample_1\tsample_2\tsample_3\tsample_4\tsample_5
    contig_1\t10\t20\t30\t40\t50
    contig_2\t11\t21\t31\t41\t51
    contig_3\t50\t5\t8\t3\t9
    contig_zero\t0\t0\t0\t0\t0
""")

ANNOTATIONS_TSV = textwrap.dedent("""\
    Contig_ID\tKO_terms\tMetaCyc_terms
    contig_1\tK00001,K00360\t
    contig_2\tK00001\t
    contig_3\tK02703\t
""")


@pytest.fixture()
def coverage_tsv(tmp_path: Path) -> Path:
    p = tmp_path / "coverage.tsv"
    p.write_text(COVERAGE_TSV)
    return p


@pytest.fixture()
def annotations_tsv(tmp_path: Path) -> Path:
    p = tmp_path / "annotations.tsv"
    p.write_text(ANNOTATIONS_TSV)
    return p


# ---------------------------------------------------------------------------
# CoAbundanceCorrelator tests
# ---------------------------------------------------------------------------

def test_correlator_finds_correlated_pair(coverage_tsv: Path):
    corr = CoAbundanceCorrelator(r_threshold=0.85, p_threshold=0.01)
    pairs = corr.correlate(coverage_tsv)
    contig_pairs = {(p["contig_a"], p["contig_b"]) for p in pairs}
    # contig_1 & contig_2 are near-perfectly correlated
    assert ("contig_1", "contig_2") in contig_pairs


def test_correlator_removes_zero_coverage(coverage_tsv: Path):
    corr = CoAbundanceCorrelator()
    pairs = corr.correlate(coverage_tsv)
    names_in_pairs = set()
    for p in pairs:
        names_in_pairs.add(p["contig_a"])
        names_in_pairs.add(p["contig_b"])
    assert "contig_zero" not in names_in_pairs


def test_correlator_missing_file():
    corr = CoAbundanceCorrelator()
    with pytest.raises(FileNotFoundError):
        corr.correlate("/no/such/coverage.tsv")


def test_correlator_too_few_samples(tmp_path: Path):
    p = tmp_path / "small.tsv"
    p.write_text("Contig_ID\tsample_1\tsample_2\ncontig_1\t10\t20\ncontig_2\t11\t21\n")
    corr = CoAbundanceCorrelator()
    with pytest.raises(ValueError, match="At least 3 samples"):
        corr.correlate(p)


def test_correlator_output_fields(coverage_tsv: Path):
    corr = CoAbundanceCorrelator(r_threshold=0.85, p_threshold=0.01)
    pairs = corr.correlate(coverage_tsv)
    assert pairs  # at least one pair
    for p in pairs:
        assert "contig_a" in p
        assert "contig_b" in p
        assert "spearman_r" in p
        assert "p_value" in p
        assert abs(p["spearman_r"]) >= 0.85
        assert p["p_value"] < 0.01


# ---------------------------------------------------------------------------
# PathwayComplementarityChecker tests
# ---------------------------------------------------------------------------

def test_pathway_checker_loads_annotations(annotations_tsv: Path):
    checker = PathwayComplementarityChecker()
    checker.load_annotations(annotations_tsv)
    # contig_1 and contig_2 both have K00001 → carbon_metabolism
    assert checker.check_pair("contig_1", "contig_2") is True
    # contig_1 (carbon_metabolism, nitrogen_metabolism) vs contig_3 (photosynthesis)
    assert checker.check_pair("contig_1", "contig_3") is False


def test_pathway_checker_missing_file():
    checker = PathwayComplementarityChecker()
    with pytest.raises(FileNotFoundError):
        checker.load_annotations("/no/such/file.tsv")


def test_pathway_checker_annotates_pairs(annotations_tsv: Path):
    checker = PathwayComplementarityChecker()
    checker.load_annotations(annotations_tsv)
    pairs = [
        {"contig_a": "contig_1", "contig_b": "contig_2", "spearman_r": 0.99, "p_value": 0.001},
        {"contig_a": "contig_1", "contig_b": "contig_3", "spearman_r": 0.90, "p_value": 0.005},
    ]
    annotated = checker.annotate_pairs(pairs)
    assert annotated[0]["metabolic_match"] is True
    assert annotated[1]["metabolic_match"] is False


# ---------------------------------------------------------------------------
# EcologicalGraphInjector tests
# ---------------------------------------------------------------------------

def test_injector_adds_edges():
    graph = nx.MultiGraph()
    graph.add_node("contig_1", length=500, node_type="unknown")
    graph.add_node("contig_2", length=800, node_type="unknown")

    injector = EcologicalGraphInjector(graph)
    pairs = [
        {
            "contig_a": "contig_1",
            "contig_b": "contig_2",
            "spearman_r": 0.91,
            "p_value": 0.003,
            "metabolic_match": True,
        }
    ]
    injector.inject(pairs)

    edges = list(graph.edges(data=True))
    assert len(edges) == 1
    _, _, attrs = edges[0]
    assert attrs["type"] == "co_abundance_guild"
    assert attrs["metabolic_match"] is True
    assert abs(attrs["weight"] - 0.91) < 1e-9


def test_injector_creates_missing_nodes():
    graph = nx.MultiGraph()
    injector = EcologicalGraphInjector(graph)
    pairs = [
        {
            "contig_a": "new_contig_1",
            "contig_b": "new_contig_2",
            "spearman_r": 0.88,
            "p_value": 0.008,
            "metabolic_match": False,
        }
    ]
    injector.inject(pairs)
    assert graph.has_node("new_contig_1")
    assert graph.has_node("new_contig_2")


# ---------------------------------------------------------------------------
# EcologicalMiner end-to-end
# ---------------------------------------------------------------------------

def test_ecological_miner_end_to_end(coverage_tsv: Path, annotations_tsv: Path):
    graph = nx.MultiGraph()
    miner = EcologicalMiner(graph=graph, r_threshold=0.85, p_threshold=0.01)
    result = miner.run(coverage_tsv=coverage_tsv, annotations_tsv=annotations_tsv)

    co_abundance_edges = [
        d for _, _, d in result.edges(data=True)
        if d.get("type") == "co_abundance_guild"
    ]
    assert len(co_abundance_edges) >= 1
    # At least one edge should connect contig_1 & contig_2
    edge_nodes = {
        (u, v)
        for u, v, d in result.edges(data=True)
        if d.get("type") == "co_abundance_guild"
    }
    assert ("contig_1", "contig_2") in edge_nodes or ("contig_2", "contig_1") in edge_nodes


def test_ecological_miner_without_annotations(coverage_tsv: Path):
    graph = nx.MultiGraph()
    miner = EcologicalMiner(graph=graph)
    result = miner.run(coverage_tsv=coverage_tsv, annotations_tsv=None)
    for _, _, d in result.edges(data=True):
        if d.get("type") == "co_abundance_guild":
            assert d["metabolic_match"] is False
