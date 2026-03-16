from pathlib import Path

from contigweaver.modules.contig_reconciler import ContigGraphReconciler
from contigweaver.modules.gfa_parser import parse_gfa


def test_reconciler_adds_segment_membership_edges_for_exact_matches(
    mini_gfa_file,
    mini_contigs_fasta,
    mini_viral_fasta,
):
    graph = parse_gfa(mini_gfa_file)
    graph.add_node("NODE_7_length_110_cov_63.7593", node_type="unknown", length=110)
    graph.add_node("NODE_37_length_124_cov_14.9354", node_type="viral", length=124)
    graph.add_edge(
        "NODE_7_length_110_cov_63.7593",
        "NODE_37_length_124_cov_14.9354",
        type="crispr_targeting",
        identity=99.0,
        coverage=1.0,
    )

    reconciler = ContigGraphReconciler(graph)
    summary = reconciler.reconcile(mini_gfa_file, [mini_contigs_fasta, mini_viral_fasta])

    assert summary["resolved_contigs"] == 2
    edge_types = [data["type"] for _, _, data in graph.edges(data=True)]
    assert edge_types.count("segment_membership") >= 2
    assert graph.has_edge("NODE_7_length_110_cov_63.7593", "7")
    assert graph.has_edge("NODE_37_length_124_cov_14.9354", "37")


def test_reconciler_supports_substring_matches(tmp_path: Path):
    gfa_path = tmp_path / "substring.gfa"
    gfa_path.write_text(
        "H\tVN:Z:1.2\n"
        "S\t42\tAAAACCCCGGGGTTTT\tDP:f:2.0\tKC:i:10\n"
        "S\t43\tTTTTGGGGCCCCAAAA\tDP:f:3.0\tKC:i:11\n"
        "L\t42\t+\t43\t+\t55M\n"
    )
    fasta_path = tmp_path / "contigs.fa"
    fasta_path.write_text(
        ">NODE_42_length_24_cov_2.0\n"
        "GGGGAAAACCCCGGGGTTTTCCCC\n"
    )

    graph = parse_gfa(gfa_path)
    graph.add_node("NODE_42_length_24_cov_2.0", node_type="unknown", length=24)
    graph.add_edge(
        "NODE_42_length_24_cov_2.0",
        "42",
        type="crispr_targeting",
        identity=99.0,
        coverage=1.0,
    )

    reconciler = ContigGraphReconciler(graph, kmer_size=8)
    summary = reconciler.reconcile(gfa_path, [fasta_path])

    assert summary["resolved_contigs"] == 1
    membership_edges = [
        data
        for u, v, data in graph.edges(data=True)
        if {u, v} == {"NODE_42_length_24_cov_2.0", "42"}
        and data.get("type") == "segment_membership"
    ]
    assert membership_edges
    assert membership_edges[0]["match_type"] == "substring"
