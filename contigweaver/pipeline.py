"""
ContigWeaver — Main Pipeline Runner
=====================================
CLI entry point that orchestrates all six modules across Stage 1 and Stage 2.

Usage (Stage 1 only):
    python -m contigweaver \\
        --gfa assembly_graph.gfa \\
        --contigs contigs.fasta \\
        --viral-contigs viral_contigs.fasta \\
        --output-dir results/

Usage (Stage 1 + Stage 2):
    python -m contigweaver \\
        --gfa assembly_graph.gfa \\
        --contigs contigs.fasta \\
        --viral-contigs viral_contigs.fasta \\
        --coverage coverage_table.tsv \\
        --annotations functional_annotations.tsv \\
        --output-dir results/
"""

import argparse
from datetime import datetime
import logging
import shlex
import sys
import uuid
from pathlib import Path
from typing import Optional

import networkx as nx

from contigweaver.modules.gfa_parser import GFAParser
from contigweaver.modules.crispr_miner import CRISPRPhageMiner
from contigweaver.modules.annotation_converter import prepare_annotations_input
from contigweaver.modules.annotation_miner import AnnotationMiner
from contigweaver.modules.binning_miner import BinningMiner
from contigweaver.modules.contig_reconciler import ContigGraphReconciler
from contigweaver.modules.graph_exporter import export_graph, export_index_report
from contigweaver.modules.ecological_miner import EcologicalMiner

logger = logging.getLogger("contigweaver")


def _detect_stage_label(args: argparse.Namespace) -> str:
    has_stage2_inputs = any(
        [
            args.coverage,
            args.annotations,
            args.annotation_data,
            args.binning,
        ]
    )
    return "stage2" if has_stage2_inputs else "stage1"


def _create_run_directory(root_dir: str | Path, stage_label: str) -> Path:
    runs_root = Path(root_dir)
    runs_root.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    unique_id = uuid.uuid4().hex[:8]
    run_dir = runs_root / f"run_{timestamp}_{stage_label}_{unique_id}"
    run_dir.mkdir(parents=True, exist_ok=False)

    latest_link = runs_root / "latest"
    tmp_link = runs_root / f".latest_tmp_{unique_id}"
    try:
        if tmp_link.exists() or tmp_link.is_symlink():
            tmp_link.unlink()
        tmp_link.symlink_to(run_dir.name, target_is_directory=True)
        tmp_link.replace(latest_link)
    except OSError as exc:
        logger.warning("Could not update latest symlink at %s: %s", latest_link, exc)
        if tmp_link.exists() or tmp_link.is_symlink():
            tmp_link.unlink()

    return run_dir


# ===========================================================================
# Pipeline class
# ===========================================================================

class ContigWeaverPipeline:
    """
    Full ContigWeaver pipeline orchestrator.

    Stage 1
    -------
    1. Module 1 — GFA Parser (physical overlap edges)
    2. Module 2 — CRISPR-Phage Miner (crispr_targeting edges)
    3. Module 3 — Graph Export (TSV + interactive HTML)

    Stage 2 (optional, requires coverage & annotations)
    -------
    4. Module 4 — Co-Abundance Correlator
    5. Module 5 — Pathway Complementarity Checker
    6. Module 6 — Graph Injection + re-export
    """

    def __init__(
        self,
        output_dir: str | Path = "contigweaver_output",
        minced_bin: str = "minced",
        blastn_bin: str = "blastn",
        makeblastdb_bin: str = "makeblastdb",
        spearman_threshold: float = 0.85,
        p_value_threshold: float = 0.01,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.minced_bin = minced_bin
        self.blastn_bin = blastn_bin
        self.makeblastdb_bin = makeblastdb_bin
        self.spearman_threshold = spearman_threshold
        self.p_value_threshold = p_value_threshold

        self.graph: nx.MultiGraph = nx.MultiGraph()
        self._stage1_gfa_path: Path | None = None
        self._stage1_fasta_paths: list[Path] = []
        self._stage2_coverage_path: Path | None = None
        self._stage2_annotations_path: Path | None = None
        self._stage2_annotation_data_path: Path | None = None
        self._stage2_binning_path: Path | None = None

    # ------------------------------------------------------------------
    # Stage 1
    # ------------------------------------------------------------------

    def run_stage1(
        self,
        gfa_path: str | Path,
        contigs_fasta: str | Path,
        viral_contigs_fasta: str | Path,
    ) -> None:
        """
        Run Stage 1: physical evidence + CRISPR-Phage evidence.

        Parameters
        ----------
        gfa_path : str | Path
            Assembly graph (.gfa).
        contigs_fasta : str | Path
            All assembled contigs (.fasta).
        viral_contigs_fasta : str | Path
            Viral contig subset (.fasta).
        """
        logger.info("=== Stage 1: Physical + CRISPR Evidence ===")
        self._stage1_gfa_path = Path(gfa_path)
        self._stage1_fasta_paths = [Path(contigs_fasta), Path(viral_contigs_fasta)]

        # Module 1 — GFA Parser
        logger.info("--- Module 1: GFA Parser ---")
        gfa_parser = GFAParser(graph=self.graph)
        gfa_parser.parse(gfa_path)
        logger.info(
            "After Module 1: %d nodes, %d edges",
            self.graph.number_of_nodes(),
            self.graph.number_of_edges(),
        )

        # Module 2 — CRISPR-Phage Miner
        logger.info("--- Module 2: CRISPR-Phage Miner ---")
        crispr_miner = CRISPRPhageMiner(
            graph=self.graph,
            minced_bin=self.minced_bin,
            blastn_bin=self.blastn_bin,
            makeblastdb_bin=self.makeblastdb_bin,
        )
        crispr_miner.run(
            contigs_fasta=contigs_fasta,
            viral_contigs_fasta=viral_contigs_fasta,
            work_dir=self.output_dir / "workdir",
        )
        logger.info(
            "After Module 2: %d nodes, %d edges",
            self.graph.number_of_nodes(),
            self.graph.number_of_edges(),
        )

        logger.info("--- Bridging contig IDs to GFA segments ---")
        self._reconcile_contig_ids()

        # Module 3 — Export
        logger.info("--- Module 3: Graph Export ---")
        viral_set = self._collect_viral_nodes()
        self._export(viral_set)

    # ------------------------------------------------------------------
    # Stage 2
    # ------------------------------------------------------------------

    def run_stage2(
        self,
        coverage_tsv: Optional[str | Path] = None,
        annotations_tsv: Optional[str | Path] = None,
        annotation_data_tsv: Optional[str | Path] = None,
        binning_tsv: Optional[str | Path] = None,
    ) -> None:
        """
        Run Stage 2: ecological co-abundance evidence (requires Stage 1 first).

        Parameters
        ----------
        coverage_tsv : str | Path
            Coverage table TSV (Contig_ID + sample columns).
        annotations_tsv : str | Path, optional
            Functional annotations TSV for pathway complementarity.
        """
        if (
            coverage_tsv is None
            and annotations_tsv is None
            and annotation_data_tsv is None
            and binning_tsv is None
        ):
            raise ValueError(
                "Stage 2 requires at least one input: --coverage, --annotations, --annotation-data, or --binning."
            )

        logger.info("=== Stage 2: Ecological + Annotation + Binning Evidence ===")
        self._stage2_coverage_path = Path(coverage_tsv) if coverage_tsv else None
        self._stage2_annotations_path = Path(annotations_tsv) if annotations_tsv else None
        self._stage2_annotation_data_path = (
            Path(annotation_data_tsv) if annotation_data_tsv else None
        )
        self._stage2_binning_path = Path(binning_tsv) if binning_tsv else None

        prepared_annotations: Optional[str | Path] = None
        if annotations_tsv is not None:
            prepared_annotations = prepare_annotations_input(
                annotations_input=annotations_tsv,
                work_dir=self.output_dir / "workdir",
            )

        if coverage_tsv is not None:
            eco_miner = EcologicalMiner(
                graph=self.graph,
                r_threshold=self.spearman_threshold,
                p_threshold=self.p_value_threshold,
            )
            eco_miner.run(
                coverage_tsv=coverage_tsv,
                annotations_tsv=prepared_annotations,
            )

        self._reconcile_contig_ids()

        annotation_input = annotation_data_tsv or prepared_annotations
        if annotation_input is not None:
            logger.info("--- Stage 2A: Annotation Miner ---")
            annotation_miner = AnnotationMiner(graph=self.graph)
            annotation_miner.run(annotation_input)

        if binning_tsv is not None:
            logger.info("--- Stage 2B: Binning Miner ---")
            binning_miner = BinningMiner(graph=self.graph)
            binning_miner.run(binning_tsv)

        logger.info(
            "After Stage 2: %d nodes, %d edges",
            self.graph.number_of_nodes(),
            self.graph.number_of_edges(),
        )

        self._reconcile_contig_ids()

        # Re-export with updated graph
        logger.info("--- Re-exporting updated graph ---")
        viral_set = self._collect_viral_nodes()
        self._export(viral_set)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _collect_viral_nodes(self) -> set[str]:
        """Return set of node IDs marked as viral in the graph."""
        return {
            n for n, d in self.graph.nodes(data=True)
            if d.get("node_type") == "viral"
        }

    def _reconcile_contig_ids(self) -> None:
        if self._stage1_gfa_path is None or not self._stage1_fasta_paths:
            return

        reconciler = ContigGraphReconciler(graph=self.graph)
        summary = reconciler.reconcile(
            gfa_path=self._stage1_gfa_path,
            fasta_paths=self._stage1_fasta_paths,
        )
        logger.info(
            "Contig bridging summary: %d candidate contigs, %d resolved, %d bridge edges.",
            summary["candidate_contigs"],
            summary["resolved_contigs"],
            summary["segment_membership_edges"],
        )

    def _export(self, viral_set: set[str]) -> None:
        tsv_path = self.output_dir / "contigweaver_edges.tsv"
        html_path = self.output_dir / "contigweaver_network.html"
        index_path = self.output_dir / "index.html"
        export_graph(
            graph=self.graph,
            tsv_path=tsv_path,
            html_path=html_path,
            viral_contigs=viral_set,
        )
        export_index_report(
            graph=self.graph,
            output_path=index_path,
            network_html_path=html_path,
            edges_tsv_path=tsv_path,
            viral_contigs=viral_set,
            input_paths=self._report_inputs(),
            related_html_paths=self._discover_related_html_reports(),
        )
        logger.info("Outputs saved to %s", self.output_dir)

    def _report_inputs(self) -> dict[str, Path | None]:
        return {
            "Assembly graph": self._stage1_gfa_path,
            "Contigs FASTA": self._stage1_fasta_paths[0] if self._stage1_fasta_paths else None,
            "Viral contigs FASTA": self._stage1_fasta_paths[1] if len(self._stage1_fasta_paths) > 1 else None,
            "Coverage TSV": self._stage2_coverage_path,
            "Annotations": self._stage2_annotations_path,
            "Annotation miner TSV": self._stage2_annotation_data_path,
            "Binning TSV": self._stage2_binning_path,
        }

    def _discover_related_html_reports(self) -> list[Path]:
        if self._stage1_gfa_path is None:
            return []

        roots = {
            self._stage1_gfa_path.parent,
            *(path.parent for path in self._stage1_fasta_paths),
        }
        related: list[Path] = []
        seen: set[Path] = set()
        for root in roots:
            for candidate_root in (root / "QC", root):
                if not candidate_root.exists():
                    continue
                for html_path in sorted(candidate_root.glob("**/*.html")):
                    if html_path in seen:
                        continue
                    seen.add(html_path)
                    related.append(html_path)
                    if len(related) >= 6:
                        return related
        return related


# ===========================================================================
# CLI
# ===========================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="contigweaver",
        description=(
            "ContigWeaver — Reconstruct inter-contig interaction networks from "
            "low-depth metagenomic data without MAG binning."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required Stage 1 inputs
    stage1 = parser.add_argument_group("Stage 1 inputs (required)")
    stage1.add_argument(
        "--gfa", required=True, metavar="FILE",
        help="Assembly graph in GFA v1 format (e.g., from MEGAHIT/metaSPAdes).",
    )
    stage1.add_argument(
        "--contigs", required=True, metavar="FILE",
        help="All assembled contigs FASTA.",
    )
    stage1.add_argument(
        "--viral-contigs", required=True, metavar="FILE",
        help="Viral contig subset FASTA (e.g., from geNomad).",
    )

    # Optional Stage 2 inputs
    stage2 = parser.add_argument_group("Stage 2 inputs (optional)")
    stage2.add_argument(
        "--coverage", metavar="FILE", default=None,
        help="Coverage table TSV (Contig_ID + per-sample coverage columns).",
    )
    stage2.add_argument(
        "--annotations", metavar="FILE", default=None,
        help=(
            "Functional annotations TSV or Prokka directory. "
            "TSV must contain Contig_ID plus KO_terms/MetaCyc_terms or functional_terms. "
            "Directories are converted to a contig-level TSV automatically for Stage 2."
        ),
    )
    stage2.add_argument(
        "--annotation-data", metavar="FILE", default=None,
        help=(
            "Annotation Miner TSV. Required column: Contig_ID. Optional columns: "
            "taxonomy_label, taxonomy_rank, taxonomy_confidence, functional_terms, "
            "cas_genes, has_crispr, spacer_count. "
            "If omitted, --annotations TSV is reused when available."
        ),
    )
    stage2.add_argument(
        "--binning", metavar="FILE", default=None,
        help="Binning TSV with Contig_ID and Bin_ID columns.",
    )

    # Output
    out_grp = parser.add_argument_group("Output")
    out_grp.add_argument(
        "--output-dir", metavar="DIR", default="runs",
        help="Root directory for per-run output folders.",
    )

    # External tool paths
    tools = parser.add_argument_group("External tool paths")
    tools.add_argument("--minced-bin", default="minced", metavar="PATH",
                       help="Path to the minced executable.")
    tools.add_argument("--blastn-bin", default="blastn", metavar="PATH",
                       help="Path to the blastn executable.")
    tools.add_argument("--makeblastdb-bin", default="makeblastdb", metavar="PATH",
                       help="Path to the makeblastdb executable.")

    # Thresholds
    thresh = parser.add_argument_group("Stage 2 thresholds")
    thresh.add_argument(
        "--spearman-threshold", type=float, default=0.85, metavar="FLOAT",
        help="Minimum |Spearman r| for co-abundance edges.",
    )
    thresh.add_argument(
        "--p-value-threshold", type=float, default=0.01, metavar="FLOAT",
        help="Maximum p-value for co-abundance edges.",
    )

    # Logging
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Enable DEBUG-level logging.",
    )

    return parser


def setup_logging(run_dir: str | Path, verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    log_path = Path(run_dir) / "run.log"

    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    root_logger.handlers.clear()

    formatter = logging.Formatter(
        fmt="%(asctime)s  %(levelname)-8s  %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(level)
    stream_handler.setFormatter(formatter)
    root_logger.addHandler(stream_handler)

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)


def main(argv: Optional[list[str]] = None) -> int:
    start_ts = datetime.now()
    parser = build_parser()
    args = parser.parse_args(argv)

    stage_label = _detect_stage_label(args)
    run_dir = _create_run_directory(args.output_dir, stage_label)
    setup_logging(run_dir=run_dir, verbose=args.verbose)

    logger.info("ContigWeaver v1.0.0 — starting pipeline.")
    logger.info("Run directory: %s", run_dir)
    logger.info("Stage mode: %s", stage_label)
    logger.info("Command: %s", " ".join(shlex.quote(arg) for arg in sys.argv))
    logger.info("Arguments: %s", vars(args))

    pipeline = ContigWeaverPipeline(
        output_dir=run_dir,
        minced_bin=args.minced_bin,
        blastn_bin=args.blastn_bin,
        makeblastdb_bin=args.makeblastdb_bin,
        spearman_threshold=args.spearman_threshold,
        p_value_threshold=args.p_value_threshold,
    )

    try:
        pipeline.run_stage1(
            gfa_path=args.gfa,
            contigs_fasta=args.contigs,
            viral_contigs_fasta=args.viral_contigs,
        )

        if args.coverage or args.annotations or args.annotation_data or args.binning:
            pipeline.run_stage2(
                coverage_tsv=args.coverage,
                annotations_tsv=args.annotations,
                annotation_data_tsv=args.annotation_data,
                binning_tsv=args.binning,
            )
        else:
            logger.info(
                "No Stage 2 inputs supplied; skipping Stage 2 (ecological/annotation/binning evidence)."
            )

    except FileNotFoundError as exc:
        logger.error("Input file error: %s", exc)
        logger.info("Pipeline failed at: %s", datetime.now().isoformat(timespec="seconds"))
        return 1
    except RuntimeError as exc:
        logger.error("Pipeline error: %s", exc)
        logger.info("Pipeline failed at: %s", datetime.now().isoformat(timespec="seconds"))
        return 1
    except Exception as exc:  # noqa: BLE001
        logger.exception("Unexpected error: %s", exc)
        logger.info("Pipeline failed at: %s", datetime.now().isoformat(timespec="seconds"))
        return 2

    end_ts = datetime.now()
    duration_s = (end_ts - start_ts).total_seconds()
    logger.info("Pipeline finished at: %s", end_ts.isoformat(timespec="seconds"))
    logger.info("Pipeline duration: %.2f seconds", duration_s)
    logger.info("Pipeline complete. Results in: %s", run_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
