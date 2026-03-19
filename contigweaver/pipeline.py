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
import csv
from datetime import datetime
import gzip
import logging
import shlex
import sys
import uuid
from pathlib import Path
from typing import Optional

import networkx as nx

from contigweaver.modules.gfa_parser import GFAParser
from contigweaver.modules.crispr_miner import CRISPRPhageMiner
from contigweaver.modules.annotation_converter import (
    METHOD_ALIASES,
    detect_method_from_text,
    detect_methods_from_inputs,
    prepare_annotations_input,
)
from contigweaver.modules.annotation_miner import AnnotationMiner
from contigweaver.modules.binning_miner import BinningMiner
from contigweaver.modules.contig_reconciler import ContigGraphReconciler
from contigweaver.modules.graph_exporter import export_graph, export_index_report
from contigweaver.modules.ecological_miner import EcologicalMiner

logger = logging.getLogger("contigweaver")
CANONICAL_BINNING_METHODS = tuple(sorted(METHOD_ALIASES.keys()))


def _parse_binning_method_selection(raw_value: str | None) -> set[str] | None:
    if raw_value is None:
        return None

    normalized = raw_value.strip().lower()
    if normalized in {"", "auto"}:
        return None
    if normalized == "all":
        return set(CANONICAL_BINNING_METHODS)

    selected: set[str] = set()
    for token in raw_value.split(","):
        item = token.strip().lower()
        if not item:
            continue
        if item == "metabat":
            item = "metabat2"
        elif item == "maxbin":
            item = "maxbin2"
        if item not in METHOD_ALIASES:
            valid = ", ".join(["auto", "all", *CANONICAL_BINNING_METHODS])
            raise ValueError(
                f"Unknown binning method '{item}' in --binning-methods. Valid values: {valid}."
            )
        selected.add(item)

    return selected or None


class CompactRepetitionFilter(logging.Filter):
    def __init__(self) -> None:
        super().__init__()
        self._last_signature: tuple[str, int, str] | None = None
        self._repeat_count = 0
        self._attached_handler: logging.Handler | None = None
        self._emitting_summary = False

    def bind_handler(self, handler: logging.Handler) -> None:
        self._attached_handler = handler

    def filter(self, record: logging.LogRecord) -> bool:
        if self._emitting_summary:
            return True

        signature = (record.name, record.levelno, record.getMessage())
        if self._last_signature is None:
            self._last_signature = signature
            return True

        if signature == self._last_signature:
            self._repeat_count += 1
            return False

        self._flush_summary(record)
        self._last_signature = signature
        return True

    def flush_pending(self) -> None:
        self._flush_summary(None)

    def _flush_summary(self, template_record: logging.LogRecord | None) -> None:
        if self._repeat_count <= 0 or self._attached_handler is None:
            self._repeat_count = 0
            return

        if template_record is None:
            summary_record = logging.LogRecord(
                name="contigweaver",
                level=logging.INFO,
                pathname=__file__,
                lineno=0,
                msg="Previous message repeated %d times.",
                args=(self._repeat_count,),
                exc_info=None,
            )
        else:
            summary_record = logging.LogRecord(
                name=template_record.name,
                level=template_record.levelno,
                pathname=template_record.pathname,
                lineno=template_record.lineno,
                msg="Previous message repeated %d times.",
                args=(self._repeat_count,),
                exc_info=None,
            )

        self._repeat_count = 0
        self._emitting_summary = True
        try:
            self._attached_handler.emit(summary_record)
        finally:
            self._emitting_summary = False


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
        self._stage2_binning_paths: list[Path] = []

    def _expand_binning_inputs(
        self,
        binning_tsv: Optional[str | Path | list[str | Path]],
    ) -> list[Path]:
        if binning_tsv is None:
            return []

        raw_inputs: list[Path]
        if isinstance(binning_tsv, (str, Path)):
            raw_inputs = [Path(binning_tsv)]
        else:
            raw_inputs = [Path(value) for value in binning_tsv]

        expanded: list[Path] = []
        for item in raw_inputs:
            if item.is_file():
                expanded.append(item)
                continue
            if item.is_dir():
                candidates = sorted(item.glob("**/*binning*.tsv"))
                detected_methods = {
                    method
                    for method in (
                        detect_method_from_text(str(path)) for path in candidates
                    )
                    if method is not None
                }
                derived = self._derive_binning_tsvs_from_bins_directory(
                    root=item,
                    skip_methods=detected_methods,
                )
                if not candidates and not derived:
                    raise FileNotFoundError(f"No binning TSV files found in directory: {item}")
                expanded.extend(candidates)
                expanded.extend(derived)
                continue
            raise FileNotFoundError(f"Binning input not found: {item}")

        deduped: list[Path] = []
        seen: set[Path] = set()
        for path in expanded:
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            deduped.append(path)
        return deduped

    @staticmethod
    def _iter_fasta_headers(path: Path) -> list[str]:
        open_func = gzip.open if str(path).endswith(".gz") else open
        headers: list[str] = []
        with open_func(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    headers.append(line[1:].strip().split()[0])
        return headers

    @staticmethod
    def _bin_id_from_fasta_path(path: Path) -> str:
        name = path.name
        for suffix in (".fasta.gz", ".fa.gz", ".fasta", ".fa"):
            if name.endswith(suffix):
                return name[: -len(suffix)]
        return path.stem

    def _derive_binning_tsvs_from_bins_directory(
        self,
        root: Path,
        skip_methods: set[str],
    ) -> list[Path]:
        method_to_fastas: dict[str, list[Path]] = {}
        for fasta_path in sorted(root.glob("**/bins/*")):
            if not fasta_path.is_file():
                continue
            if not any(str(fasta_path).endswith(ext) for ext in (".fa", ".fasta", ".fa.gz", ".fasta.gz")):
                continue
            method = detect_method_from_text(str(fasta_path))
            if method is None or method in skip_methods:
                continue
            method_to_fastas.setdefault(method, []).append(fasta_path)

        if not method_to_fastas:
            return []

        work_dir = self.output_dir / "workdir"
        work_dir.mkdir(parents=True, exist_ok=True)
        generated: list[Path] = []
        for method, fasta_files in sorted(method_to_fastas.items()):
            out_tsv = work_dir / f"{method}_auto_binning.tsv"
            with out_tsv.open("w", newline="") as fh:
                writer = csv.writer(fh, delimiter="\t")
                writer.writerow(["Contig_ID", "Bin_ID"])
                for fasta_path in fasta_files:
                    bin_id = self._bin_id_from_fasta_path(fasta_path)
                    for contig_id in self._iter_fasta_headers(fasta_path):
                        writer.writerow([contig_id, bin_id])
            logger.info(
                "Derived binning TSV for method '%s' from %d FASTA bin files: %s",
                method,
                len(fasta_files),
                out_tsv,
            )
            generated.append(out_tsv)

        return generated

    @staticmethod
    def _filter_binning_paths_by_methods(
        binning_paths: list[Path],
        selected_methods: set[str] | None,
    ) -> list[Path]:
        if not selected_methods:
            return binning_paths

        filtered: list[Path] = []
        for path in binning_paths:
            method = detect_method_from_text(str(path))
            if method is None or method in selected_methods:
                filtered.append(path)

        if filtered:
            return filtered

        raise ValueError(
            "No binning TSV files remained after --binning-methods filtering. "
            f"Requested methods: {', '.join(sorted(selected_methods))}."
        )

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
        binning_tsv: Optional[str | Path | list[str | Path]] = None,
        binning_methods: Optional[str] = None,
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
        resolved_binning_paths = self._expand_binning_inputs(binning_tsv)
        selected_methods = _parse_binning_method_selection(binning_methods)
        selected_binning_paths = self._filter_binning_paths_by_methods(
            resolved_binning_paths,
            selected_methods,
        )

        self._stage2_coverage_path = Path(coverage_tsv) if coverage_tsv else None
        self._stage2_annotations_path = Path(annotations_tsv) if annotations_tsv else None
        self._stage2_annotation_data_path = (
            Path(annotation_data_tsv) if annotation_data_tsv else None
        )
        self._stage2_binning_paths = selected_binning_paths

        annotation_method_selection = selected_methods
        if annotation_method_selection is None:
            annotation_method_selection = detect_methods_from_inputs(
                [str(path) for path in selected_binning_paths]
            )

        prepared_annotations: Optional[str | Path] = None
        if annotations_tsv is not None:
            prepared_annotations = prepare_annotations_input(
                annotations_input=annotations_tsv,
                work_dir=self.output_dir / "workdir",
                binning_input=[str(path) for path in selected_binning_paths],
                method_selection=annotation_method_selection,
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

        if selected_binning_paths:
            logger.info("--- Stage 2B: Binning Miner ---")
            binning_miner = BinningMiner(graph=self.graph)
            for binning_path in selected_binning_paths:
                logger.info("Applying binning assignments from %s", binning_path)
                binning_miner.run(binning_path)

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
        report_inputs: dict[str, Path | None] = {
            "Assembly graph": self._stage1_gfa_path,
            "Contigs FASTA": self._stage1_fasta_paths[0] if self._stage1_fasta_paths else None,
            "Viral contigs FASTA": self._stage1_fasta_paths[1] if len(self._stage1_fasta_paths) > 1 else None,
            "Coverage TSV": self._stage2_coverage_path,
            "Annotations": self._stage2_annotations_path,
            "Annotation miner TSV": self._stage2_annotation_data_path,
        }
        if len(self._stage2_binning_paths) == 1:
            report_inputs["Binning TSV"] = self._stage2_binning_paths[0]
        elif self._stage2_binning_paths:
            for index, path in enumerate(self._stage2_binning_paths, start=1):
                report_inputs[f"Binning TSV #{index}"] = path
        return report_inputs

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
            "Directories are converted to a contig-level TSV automatically for Stage 2 and "
            "matched to --binning method when method tags are detectable in names."
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
        "--binning", metavar="FILE_OR_DIR", nargs="+", default=None,
        help=(
            "One or more binning TSV files (Contig_ID, Bin_ID) or directories. "
            "Directories are scanned recursively for *binning*.tsv files."
        ),
    )
    stage2.add_argument(
        "--binning-methods", metavar="METHODS", default="auto",
        help=(
            "Method selection for binning+annotation matching: auto (default), all, or "
            f"comma-separated subset from {', '.join(CANONICAL_BINNING_METHODS)}."
        ),
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


def setup_logging(run_dir: str | Path, verbose: bool = False) -> list[CompactRepetitionFilter]:
    level = logging.DEBUG if verbose else logging.INFO
    log_path = Path(run_dir) / "run.log"

    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    root_logger.handlers.clear()

    formatter = logging.Formatter(
        fmt="%(asctime)s %(levelname).1s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(level)
    stream_handler.setFormatter(formatter)
    stream_filter = CompactRepetitionFilter()
    stream_filter.bind_handler(stream_handler)
    stream_handler.addFilter(stream_filter)
    root_logger.addHandler(stream_handler)

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    file_filter = CompactRepetitionFilter()
    file_filter.bind_handler(file_handler)
    file_handler.addFilter(file_filter)
    root_logger.addHandler(file_handler)

    for noisy_logger in ("matplotlib", "numexpr", "urllib3", "asyncio"):
        logging.getLogger(noisy_logger).setLevel(logging.WARNING)

    return [stream_filter, file_filter]


def main(argv: Optional[list[str]] = None) -> int:
    start_ts = datetime.now()
    parser = build_parser()
    args = parser.parse_args(argv)

    stage_label = _detect_stage_label(args)
    run_dir = _create_run_directory(args.output_dir, stage_label)
    repetition_filters = setup_logging(run_dir=run_dir, verbose=args.verbose)

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
                binning_methods=args.binning_methods,
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

    for repetition_filter in repetition_filters:
        repetition_filter.flush_pending()
    for handler in logging.getLogger().handlers:
        handler.flush()
    return 0


if __name__ == "__main__":
    sys.exit(main())
