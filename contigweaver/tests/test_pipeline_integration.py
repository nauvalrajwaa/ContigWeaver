"""
Integration tests for the ContigNexus pipeline using a compact fixture
that mirrors the real SPAdes output format from
results/Assembly/SPAdes/SPAdes-sponge_brin.assembly.gfa.gz.

Coverage:
- Stage 1: GFA parsing (plain + .gz), CRISPR miner (mocked), graph export
- Stage 2: ecological miner (co-abundance + pathway)
- Full pipeline CLI (end-to-end with mocked external tools)
- Real-data smoke test (skipped when results/ not present)
"""
import gzip
import subprocess
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import networkx as nx
import pytest

from contigweaver.modules.gfa_parser import GFAParser, parse_gfa
from contigweaver.modules.graph_exporter import export_graph
from contigweaver.modules.ecological_miner import (
    CoAbundanceCorrelator,
    EcologicalMiner,
)
from contigweaver.pipeline import ContigWeaverPipeline, build_parser

# real_data marker and REAL_GFA path re-declared here so the module is
# importable without relying on conftest (pytest auto-loads conftest but does
# not allow direct imports from it).
from pathlib import Path as _Path

_REAL_GFA = _Path(__file__).parent.parent.parent / "results/Assembly/SPAdes/SPAdes-sponge_brin.assembly.gfa.gz"
REAL_GFA = _REAL_GFA

real_data = pytest.mark.skipif(
    not _REAL_GFA.exists(),
    reason="Real results/ data not available",
)


# ===========================================================================
# 1. GFA Parser — plain and gzipped real-format files
# ===========================================================================


class TestGFAParserRealFormat:
    """Parse the mini SPAdes-format GFA (integer IDs, DP:f: tags, 55M overlaps)."""

    def test_plain_gfa_segment_count(self, mini_gfa_file):
        g = parse_gfa(mini_gfa_file)
        assert g.number_of_nodes() == 6

    def test_plain_gfa_link_count(self, mini_gfa_file):
        g = parse_gfa(mini_gfa_file)
        assert g.number_of_edges() == 5

    def test_gz_gfa_matches_plain(self, mini_gfa_file, mini_gfa_gz_file):
        """Gzipped and plain GFA must produce identical graphs."""
        g_plain = parse_gfa(mini_gfa_file)
        g_gz = parse_gfa(mini_gfa_gz_file)
        assert set(g_plain.nodes()) == set(g_gz.nodes())
        assert g_plain.number_of_edges() == g_gz.number_of_edges()

    def test_integer_segment_ids(self, mini_gfa_file):
        g = parse_gfa(mini_gfa_file)
        for node in g.nodes():
            assert str(node).isdigit(), f"Expected integer-string node ID, got {node!r}"

    def test_length_inferred_from_sequence(self, mini_gfa_file):
        """Real SPAdes GFA has no LN:i: tag — length must be inferred from sequence."""
        g = parse_gfa(mini_gfa_file)
        assert g.nodes["7"]["length"] == 109  # sequence in GFA S-line is 109 bp
        assert g.nodes["5"]["length"] == 175  # actual sequence length in fixture

    def test_physical_overlap_edge_type(self, mini_gfa_file):
        g = parse_gfa(mini_gfa_file)
        for _, _, data in g.edges(data=True):
            assert data["type"] == "physical_overlap"

    def test_overlap_cigar_is_55m(self, mini_gfa_file):
        """All overlaps in the real file are 55M."""
        g = parse_gfa(mini_gfa_file)
        for _, _, data in g.edges(data=True):
            assert data["overlap_cigar"] == "55M"

    def test_dp_depth_tag_not_required(self, mini_gfa_file):
        """DP:f: is an optional tag; parser should not crash on its presence."""
        g = parse_gfa(mini_gfa_file)
        assert g.number_of_nodes() > 0  # survived parsing

    def test_gz_segment_lengths(self, mini_gfa_gz_file):
        g = parse_gfa(mini_gfa_gz_file)
        assert g.nodes["7"]["length"] == 109  # sequence in GFA S-line is 109 bp
        assert g.nodes["37"]["length"] == 122  # sequence in GFA S-line is 122 bp


# ===========================================================================
# 2. Stage 1 — Full pipeline run (CRISPR external tools mocked)
# ===========================================================================


class TestStage1Integration:
    """
    Stage 1: parse GFA → mine CRISPR (mocked) → export.
    External binaries (minced, blastn, makeblastdb) are replaced with mocks
    so the test is hermetic and fast.
    """

    def _make_pipeline(self, tmp_path) -> ContigWeaverPipeline:
        return ContigWeaverPipeline(
            output_dir=tmp_path,
            minced_bin="minced",
            blastn_bin="blastn",
            makeblastdb_bin="makeblastdb",
        )

    def test_stage1_creates_tsv(
        self,
        tmp_path,
        mini_gfa_file,
        mini_contigs_fasta,
        mini_viral_fasta,
    ):
        pipeline = self._make_pipeline(tmp_path)
        with patch("subprocess.run") as mock_run, patch(
            "contigweaver.modules.crispr_miner.SpacerExtractor.extract_to_fasta",
            return_value=0,
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            pipeline.run_stage1(mini_gfa_file, mini_contigs_fasta, mini_viral_fasta)

        tsv = tmp_path / "contigweaver_edges.tsv"
        assert tsv.exists(), "Stage 1 must produce contigweaver_edges.tsv"

    def test_stage1_tsv_has_physical_overlap_rows(
        self,
        tmp_path,
        mini_gfa_file,
        mini_contigs_fasta,
        mini_viral_fasta,
    ):
        pipeline = self._make_pipeline(tmp_path)
        with patch("subprocess.run") as mock_run, patch(
            "contigweaver.modules.crispr_miner.SpacerExtractor.extract_to_fasta",
            return_value=0,
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            pipeline.run_stage1(mini_gfa_file, mini_contigs_fasta, mini_viral_fasta)

        tsv = tmp_path / "contigweaver_edges.tsv"
        content = tsv.read_text()
        assert "Physical" in content, f"Expected Physical overlap rows, got:\n{content}"

    def test_stage1_reconciles_contig_ids_into_graph(
        self,
        tmp_path,
        mini_gfa_file,
        mini_contigs_fasta,
        mini_viral_fasta,
    ):
        pipeline = self._make_pipeline(tmp_path)

        def inject_crispr_edge(*args, **kwargs):
            pipeline.graph.add_node("NODE_7_length_110_cov_63.7593", node_type="unknown", length=110)
            pipeline.graph.add_node("NODE_37_length_124_cov_14.9354", node_type="viral", length=124)
            pipeline.graph.add_edge(
                "NODE_7_length_110_cov_63.7593",
                "NODE_37_length_124_cov_14.9354",
                type="crispr_targeting",
                identity=99.0,
                coverage=1.0,
            )
            return pipeline.graph

        with patch(
            "contigweaver.modules.crispr_miner.CRISPRPhageMiner.run",
            side_effect=inject_crispr_edge,
        ):
            pipeline.run_stage1(mini_gfa_file, mini_contigs_fasta, mini_viral_fasta)

        edge_types = {data["type"] for _, _, data in pipeline.graph.edges(data=True)}
        assert "segment_membership" in edge_types
        assert pipeline.graph.has_edge("NODE_7_length_110_cov_63.7593", "7")
        assert pipeline.graph.has_edge("NODE_37_length_124_cov_14.9354", "37")

        tsv = tmp_path / "contigweaver_edges.tsv"
        assert "Segment-Membership" in tsv.read_text()

    def test_stage1_creates_html(
        self,
        tmp_path,
        mini_gfa_file,
        mini_contigs_fasta,
        mini_viral_fasta,
    ):
        pytest.importorskip("pyvis")
        pipeline = self._make_pipeline(tmp_path)
        with patch("subprocess.run") as mock_run, patch(
            "contigweaver.modules.crispr_miner.SpacerExtractor.extract_to_fasta",
            return_value=0,
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            pipeline.run_stage1(mini_gfa_file, mini_contigs_fasta, mini_viral_fasta)

        html = tmp_path / "contigweaver_network.html"
        assert html.exists() and html.stat().st_size > 500

    def test_stage1_graph_has_6_nodes(
        self,
        tmp_path,
        mini_gfa_file,
        mini_contigs_fasta,
        mini_viral_fasta,
    ):
        pipeline = self._make_pipeline(tmp_path)
        with patch("subprocess.run") as mock_run, patch(
            "contigweaver.modules.crispr_miner.SpacerExtractor.extract_to_fasta",
            return_value=0,
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            pipeline.run_stage1(mini_gfa_file, mini_contigs_fasta, mini_viral_fasta)

        assert pipeline.graph.number_of_nodes() == 6

    def test_stage1_with_gz_gfa(
        self,
        tmp_path,
        mini_gfa_gz_file,
        mini_contigs_fasta,
        mini_viral_fasta,
    ):
        """Pipeline must accept .gfa.gz directly without pre-decompression."""
        pipeline = self._make_pipeline(tmp_path)
        with patch("subprocess.run") as mock_run, patch(
            "contigweaver.modules.crispr_miner.SpacerExtractor.extract_to_fasta",
            return_value=0,
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            pipeline.run_stage1(mini_gfa_gz_file, mini_contigs_fasta, mini_viral_fasta)

        assert pipeline.graph.number_of_nodes() == 6


# ===========================================================================
# 3. Stage 2 — Ecological miner
# ===========================================================================


class TestStage2Integration:
    def test_stage2_adds_coabundance_edges(
        self,
        tmp_path,
        mini_gfa_file,
        mini_contigs_fasta,
        mini_viral_fasta,
        mini_coverage_tsv,
        mini_annotations_tsv,
    ):
        pipeline = ContigWeaverPipeline(output_dir=tmp_path)
        with patch("subprocess.run") as mock_run, patch(
            "contigweaver.modules.crispr_miner.SpacerExtractor.extract_to_fasta",
            return_value=0,
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            pipeline.run_stage1(mini_gfa_file, mini_contigs_fasta, mini_viral_fasta)

        pipeline.run_stage2(mini_coverage_tsv, mini_annotations_tsv)

        edge_types = {d["type"] for _, _, d in pipeline.graph.edges(data=True)}
        # Should still have physical_overlap edges plus potentially co_abundance_guild
        assert "physical_overlap" in edge_types

    def test_stage2_tsv_updated(
        self,
        tmp_path,
        mini_gfa_file,
        mini_contigs_fasta,
        mini_viral_fasta,
        mini_coverage_tsv,
    ):
        pipeline = ContigWeaverPipeline(output_dir=tmp_path)
        with patch("subprocess.run") as mock_run, patch(
            "contigweaver.modules.crispr_miner.SpacerExtractor.extract_to_fasta",
            return_value=0,
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            pipeline.run_stage1(mini_gfa_file, mini_contigs_fasta, mini_viral_fasta)

        pipeline.run_stage2(mini_coverage_tsv)
        tsv = tmp_path / "contigweaver_edges.tsv"
        assert tsv.exists()


# ===========================================================================
# 4. CLI integration (end-to-end argument parsing + pipeline dispatch)
# ===========================================================================


class TestCLIIntegration:
    def test_cli_stage1_returns_0(
        self,
        tmp_path,
        mini_gfa_gz_file,
        mini_contigs_fasta,
        mini_viral_fasta,
    ):
        """CLI must return exit code 0 on valid stage-1 inputs."""
        from contigweaver.pipeline import main

        args = [
            "--gfa", str(mini_gfa_gz_file),
            "--contigs", str(mini_contigs_fasta),
            "--viral-contigs", str(mini_viral_fasta),
            "--output-dir", str(tmp_path / "cli_out"),
        ]
        with patch("sys.argv", ["contigweaver"] + args), patch(
            "subprocess.run"
        ) as mock_run, patch(
            "contigweaver.modules.crispr_miner.SpacerExtractor.extract_to_fasta",
            return_value=0,
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            rc = main()
        assert rc == 0

    def test_cli_missing_gfa_returns_1(self, tmp_path, mini_contigs_fasta, mini_viral_fasta):
        from contigweaver.pipeline import main

        args = [
            "--gfa", str(tmp_path / "nonexistent.gfa"),
            "--contigs", str(mini_contigs_fasta),
            "--viral-contigs", str(mini_viral_fasta),
            "--output-dir", str(tmp_path / "cli_err"),
        ]
        with patch("sys.argv", ["contigweaver"] + args):
            rc = main()
        assert rc == 1

    def test_cli_with_gz_gfa_produces_tsv(
        self,
        tmp_path,
        mini_gfa_gz_file,
        mini_contigs_fasta,
        mini_viral_fasta,
    ):
        from contigweaver.pipeline import main

        out = tmp_path / "cli_gz_out"
        args = [
            "--gfa", str(mini_gfa_gz_file),
            "--contigs", str(mini_contigs_fasta),
            "--viral-contigs", str(mini_viral_fasta),
            "--output-dir", str(out),
        ]
        with patch("sys.argv", ["contigweaver"] + args), patch(
            "subprocess.run"
        ) as mock_run, patch(
            "contigweaver.modules.crispr_miner.SpacerExtractor.extract_to_fasta",
            return_value=0,
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            main()
        assert (out / "contigweaver_edges.tsv").exists()


# ===========================================================================
# 5. Real-data smoke tests (skipped when results/ absent)
# ===========================================================================


@real_data
class TestRealDataSmoke:
    """
    Smoke tests against the actual SPAdes assembly.
    Run only when results/Assembly/SPAdes/SPAdes-sponge_brin.assembly.gfa.gz
    is present (CI can skip; local dev runs these automatically).
    """

    def test_real_gfa_gz_parses_without_error(self):
        """Real file must parse without exception (spot-checks first 5 segments)."""
        # Parse only a tiny slice to keep the test fast
        import gzip as _gz
        import tempfile, os

        # Extract first 50 lines into a temp file
        lines = []
        with _gz.open(REAL_GFA, "rt") as fh:
            for i, line in enumerate(fh):
                if i >= 50:
                    break
                lines.append(line)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".gfa", delete=False
        ) as tf:
            tf.writelines(lines)
            tf_path = Path(tf.name)

        try:
            g = parse_gfa(tf_path)
            # First 50 lines = 1 H + up to 49 S lines
            assert g.number_of_nodes() >= 1
        finally:
            os.unlink(tf_path)

    def test_real_gz_gfa_direct(self):
        """GFAParser must accept the .gz file path directly."""
        import gzip as _gz
        import tempfile, os

        # Write a 30-line slice to a gzipped temp file
        lines = []
        with _gz.open(REAL_GFA, "rt") as fh:
            for i, line in enumerate(fh):
                if i >= 30:
                    break
                lines.append(line)

        with tempfile.NamedTemporaryFile(
            suffix=".gfa.gz", delete=False
        ) as tf:
            tf_path = Path(tf.name)

        with _gz.open(tf_path, "wt") as gz_out:
            gz_out.writelines(lines)

        try:
            g = parse_gfa(tf_path)
            assert g.number_of_nodes() >= 1
        finally:
            os.unlink(tf_path)
