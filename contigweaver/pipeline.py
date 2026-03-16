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
import logging
import sys
from pathlib import Path
from typing import Optional

import networkx as nx

from contigweaver.modules.gfa_parser import GFAParser
from contigweaver.modules.crispr_miner import CRISPRPhageMiner
from contigweaver.modules.contig_reconciler import ContigGraphReconciler
from contigweaver.modules.graph_exporter import export_graph
from contigweaver.modules.ecological_miner import EcologicalMiner

logger = logging.getLogger("contigweaver")


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
        coverage_tsv: str | Path,
        annotations_tsv: Optional[str | Path] = None,
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
        logger.info("=== Stage 2: Ecological Evidence ===")

        eco_miner = EcologicalMiner(
            graph=self.graph,
            r_threshold=self.spearman_threshold,
            p_threshold=self.p_value_threshold,
        )
        eco_miner.run(
            coverage_tsv=coverage_tsv,
            annotations_tsv=annotations_tsv,
        )

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
        export_graph(
            graph=self.graph,
            tsv_path=tsv_path,
            html_path=html_path,
            viral_contigs=viral_set,
        )
        logger.info("Outputs saved to %s", self.output_dir)


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
            "Functional annotations TSV (Contig_ID, KO_terms, MetaCyc_terms). "
            "Required for pathway complementarity check in Stage 2."
        ),
    )

    # Output
    out_grp = parser.add_argument_group("Output")
    out_grp.add_argument(
        "--output-dir", metavar="DIR", default="contigweaver_output",
        help="Directory for all outputs.",
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


def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s  %(levelname)-8s  %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)

    logger.info("ContigWeaver v1.0.0 — starting pipeline.")

    pipeline = ContigWeaverPipeline(
        output_dir=args.output_dir,
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

        if args.coverage:
            pipeline.run_stage2(
                coverage_tsv=args.coverage,
                annotations_tsv=args.annotations,
            )
        else:
            logger.info(
                "No --coverage supplied; skipping Stage 2 (ecological evidence)."
            )

    except FileNotFoundError as exc:
        logger.error("Input file error: %s", exc)
        return 1
    except RuntimeError as exc:
        logger.error("Pipeline error: %s", exc)
        return 1
    except Exception as exc:  # noqa: BLE001
        logger.exception("Unexpected error: %s", exc)
        return 2

    logger.info("Pipeline complete. Results in: %s", args.output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
