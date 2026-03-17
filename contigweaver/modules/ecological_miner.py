"""
Module 4 & 5 & 6: Ecological Evidence Layer (Stage 2)
======================================================
Module 4 — Co-Abundance Correlator
    Reads a coverage table (contigs × samples), computes pairwise Spearman
    correlations, and retains pairs with r >= 0.85 and p < 0.01.

Module 5 — Pathway Complementarity Checker
    For each correlated pair, checks functional_annotations.tsv for shared
    macro-pathway membership using KEGG KO / MetaCyc terms.

Module 6 — Graph Injection
    Adds co_abundance_guild edges (with metabolic_match flag) to the
    existing NetworkX graph and re-exports outputs.

The three modules are encapsulated in the ``EcologicalMiner`` class.
"""

import logging
from itertools import combinations  # noqa: F401 (kept for optional use)
from pathlib import Path
from typing import Optional

import networkx as nx
import pandas as pd
from scipy.stats import spearmanr  # type: ignore

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SPEARMAN_THRESHOLD = 0.85   # minimum |r|
P_VALUE_THRESHOLD = 0.01    # maximum p-value

# Macro-pathway groupings: KO prefix / MetaCyc term → macro-pathway label.
# Extend this dictionary as needed.
MACRO_PATHWAY_MAP: dict[str, str] = {
    # Carbon metabolism
    "K00001": "carbon_metabolism", "K00002": "carbon_metabolism",
    "K00003": "carbon_metabolism", "K00016": "carbon_metabolism",
    "K00024": "carbon_metabolism", "K00025": "carbon_metabolism",
    "K00026": "carbon_metabolism", "K00036": "carbon_metabolism",
    "K00058": "carbon_metabolism", "K00615": "carbon_metabolism",
    # Nitrogen metabolism
    "K00360": "nitrogen_metabolism", "K00362": "nitrogen_metabolism",
    "K00363": "nitrogen_metabolism", "K01915": "nitrogen_metabolism",
    "K01923": "nitrogen_metabolism", "K02586": "nitrogen_metabolism",
    "K02591": "nitrogen_metabolism", "K10534": "nitrogen_metabolism",
    # Sulfur metabolism
    "K00392": "sulfur_metabolism", "K00394": "sulfur_metabolism",
    "K00395": "sulfur_metabolism", "K11180": "sulfur_metabolism",
    "K11181": "sulfur_metabolism",
    # Phosphorus metabolism
    "K01077": "phosphorus_metabolism", "K03307": "phosphorus_metabolism",
    # Photosynthesis
    "K02703": "photosynthesis", "K02704": "photosynthesis",
    "K02705": "photosynthesis",
    # Fatty acid metabolism
    "K00059": "fatty_acid_metabolism", "K00209": "fatty_acid_metabolism",
    "K01961": "fatty_acid_metabolism",
    # MetaCyc pathway IDs (partial list)
    "PWY-5505": "carbon_metabolism", "GLYCOLYSIS": "carbon_metabolism",
    "TCA": "carbon_metabolism", "NITRATE-DEG": "nitrogen_metabolism",
    "SULFATE-REDUCTION": "sulfur_metabolism", "PHOTOSYSTEM": "photosynthesis",
}


# ===========================================================================
# Module 4 — Co-Abundance Correlator
# ===========================================================================

class CoAbundanceCorrelator:
    """
    Compute pairwise Spearman correlations from a contig coverage table.

    The coverage table is a TSV with columns:
        Contig_ID  |  sample_1  |  sample_2  |  …  |  sample_N

    Zero-coverage contigs (all samples = 0) are removed before correlation.
    """

    def __init__(
        self,
        r_threshold: float = SPEARMAN_THRESHOLD,
        p_threshold: float = P_VALUE_THRESHOLD,
    ):
        self.r_threshold = r_threshold
        self.p_threshold = p_threshold

    def correlate(self, coverage_tsv: str | Path) -> list[dict]:
        """
        Read *coverage_tsv* and return correlated pairs.

        Parameters
        ----------
        coverage_tsv : str | Path
            Path to the coverage table TSV.

        Returns
        -------
        list[dict]
            Each dict: {contig_a, contig_b, spearman_r, p_value}

        Raises
        ------
        FileNotFoundError
            If *coverage_tsv* does not exist.
        ValueError
            If the table has fewer than 3 samples (not enough for correlation).
        """
        coverage_tsv = Path(coverage_tsv)
        if not coverage_tsv.exists():
            raise FileNotFoundError(
                f"Coverage table not found: {coverage_tsv}"
            )

        df = pd.read_csv(coverage_tsv, sep="\t", index_col=0)

        n_samples = df.shape[1]
        if n_samples < 3:
            raise ValueError(
                f"Coverage table has only {n_samples} sample column(s). "
                "At least 3 samples are required for meaningful Spearman correlation."
            )

        # Remove zero-coverage contigs
        non_zero_mask = (df > 0).any(axis=1)
        n_removed = (~non_zero_mask).sum()
        if n_removed > 0:
            logger.info(
                "Removed %d zero-coverage contigs before correlation.", n_removed
            )
        df = df[non_zero_mask]

        contigs = df.index.tolist()
        logger.info(
            "Computing Spearman correlations for %d contigs × %d samples …",
            len(contigs), n_samples,
        )

        pairs: list[dict] = []

        # ------------------------------------------------------------------
        # Vectorized pairwise Spearman via scipy bulk call.
        # spearmanr(matrix, axis=1) computes all n*(n-1)/2 pairs at once,
        # returning an (n, n) correlation matrix and matching p-value matrix.
        # This is dramatically faster than the O(n²) combinations loop for
        # large contig counts (tens of thousands).
        # ------------------------------------------------------------------
        import numpy as np
        from scipy.stats import spearmanr  # type: ignore

        values = df.values  # shape (n_contigs, n_samples)
        n = len(contigs)

        # scipy.stats.spearmanr on a 2-D array with axis=1 treats each row
        # as an observation vector and returns (correlation_matrix, p_matrix).
        # For n_contigs rows and n_samples columns we call it on the
        # *transposed* array so each row is a sample.
        corr_result = spearmanr(values, axis=1)  # correlate rows (contigs)
        corr_matrix = np.atleast_2d(corr_result.statistic) if hasattr(corr_result, 'statistic') else np.atleast_2d(corr_result.correlation)
        pval_matrix = np.atleast_2d(corr_result.pvalue)

        # Iterate only upper triangle (i < j) to avoid duplicates
        rows_i, cols_j = np.triu_indices(n, k=1)
        r_vals = corr_matrix[rows_i, cols_j]
        p_vals = pval_matrix[rows_i, cols_j]

        mask = (np.abs(r_vals) >= self.r_threshold) & (p_vals < self.p_threshold)
        for idx in np.nonzero(mask)[0]:
            i, j = int(rows_i[idx]), int(cols_j[idx])
            pairs.append(
                {
                    "contig_a": contigs[i],
                    "contig_b": contigs[j],
                    "spearman_r": float(r_vals[idx]),
                    "p_value": float(p_vals[idx]),
                }
            )

        logger.info(
            "Found %d correlated contig pairs (|r| >= %.2f, p < %.3f).",
            len(pairs), self.r_threshold, self.p_threshold,
        )
        return pairs


# ===========================================================================
# Module 5 — Pathway Complementarity Checker
# ===========================================================================

class PathwayComplementarityChecker:
    """
    Check whether two contigs share a macro-pathway based on their
    functional annotations (KEGG KO / MetaCyc terms).

    Annotations table TSV format (two accepted layouts):

    **Two-column layout** (KO and MetaCyc in separate columns):
        Contig_ID  |  KO_terms  |  MetaCyc_terms

    **Single-column layout** (spec default — mixed KO/MetaCyc):
        Contig_ID  |  functional_terms

    In both cases the term values are comma-separated lists.
    All columns other than Contig_ID are treated as term columns when
    neither 'KO_terms' nor 'MetaCyc_terms' headers are present.
    """
    def __init__(
        self,
        pathway_map: Optional[dict[str, str]] = None,
    ):
        self._map = pathway_map if pathway_map is not None else MACRO_PATHWAY_MAP
        self._contig_pathways: dict[str, set[str]] = {}

    def load_annotations(self, annotations_tsv: str | Path) -> None:
        """
        Load the functional annotations table.

        Parameters
        ----------
        annotations_tsv : str | Path
            Path to the functional annotations TSV.

        Raises
        ------
        FileNotFoundError
            If the file is missing.
        """
        annotations_tsv = Path(annotations_tsv)
        if not annotations_tsv.exists():
            raise FileNotFoundError(
                f"Functional annotations file not found: {annotations_tsv}"
            )

        df = pd.read_csv(annotations_tsv, sep="\t")
        # Expected columns: Contig_ID, KO_terms, MetaCyc_terms (others ignored)
        required_cols = {"Contig_ID"}
        missing = required_cols - set(df.columns)
        if missing:
            raise ValueError(
                f"Annotations table missing required columns: {missing}\n"
                f"Found: {list(df.columns)}"
            )

        for _, row in df.iterrows():
            contig_id = str(row["Contig_ID"])
            pathways: set[str] = set()

            # Support both two-column layout (KO_terms / MetaCyc_terms) and
            # single-column layout (any other non-ID column, e.g. 'functional_terms').
            known_term_cols = {"KO_terms", "MetaCyc_terms"}
            term_columns = [
                c for c in df.columns
                if c != "Contig_ID" and (
                    c in known_term_cols
                    or not known_term_cols.intersection(df.columns)  # single-col layout
                )
            ]

            for col in term_columns:
                if pd.notna(row[col]):
                    terms = [t.strip() for t in str(row[col]).split(",")]
                    for term in terms:
                        pathway = self._map.get(term)
                        if pathway:
                            pathways.add(pathway)

            self._contig_pathways[contig_id] = pathways

        logger.info(
            "Loaded pathway annotations for %d contigs.", len(self._contig_pathways)
        )

    def check_pair(self, contig_a: str, contig_b: str) -> bool:
        """
        Return True if *contig_a* and *contig_b* share at least one macro-pathway.
        """
        pathways_a = self._contig_pathways.get(contig_a, set())
        pathways_b = self._contig_pathways.get(contig_b, set())
        return bool(pathways_a & pathways_b)

    def annotate_pairs(self, pairs: list[dict]) -> list[dict]:
        """
        Add ``metabolic_match`` boolean to each pair in *pairs*.

        Parameters
        ----------
        pairs : list[dict]
            Output from :class:`CoAbundanceCorrelator` (must have
            ``contig_a`` and ``contig_b`` keys).

        Returns
        -------
        list[dict]
            Same list, each dict extended with ``metabolic_match``.
        """
        for pair in pairs:
            pair["metabolic_match"] = self.check_pair(
                pair["contig_a"], pair["contig_b"]
            )
        metabolic_matches = sum(1 for p in pairs if p["metabolic_match"])
        logger.info(
            "Pathway complementarity: %d / %d pairs have metabolic_match=True.",
            metabolic_matches, len(pairs),
        )
        return pairs


# ===========================================================================
# Module 6 — Graph Injection (Stage 2)
# ===========================================================================

class EcologicalGraphInjector:
    """Add co_abundance_guild edges from annotated pairs into a NetworkX graph."""

    def __init__(self, graph: nx.MultiGraph):
        self.graph = graph

    def inject(self, annotated_pairs: list[dict]) -> nx.MultiGraph:
        """
        Add *co_abundance_guild* edges to ``self.graph``.

        Parameters
        ----------
        annotated_pairs : list[dict]
            Dicts with keys: contig_a, contig_b, spearman_r, p_value,
            metabolic_match.

        Returns
        -------
        nx.MultiGraph
            Updated graph.
        """
        added = 0
        for pair in annotated_pairs:
            a, b = pair["contig_a"], pair["contig_b"]

            for node in (a, b):
                if not self.graph.has_node(node):
                    self.graph.add_node(node, length=0, node_type="unknown")

            self.graph.add_edge(
                a, b,
                type="co_abundance_guild",
                weight=pair["spearman_r"],
                p_value=pair["p_value"],
                metabolic_match=pair["metabolic_match"],
            )
            added += 1

        logger.info("Injected %d co_abundance_guild edges into graph.", added)
        return self.graph


# ===========================================================================
# High-Level EcologicalMiner (all three modules)
# ===========================================================================

class EcologicalMiner:
    """
    Orchestrates Modules 4, 5, and 6 for Stage 2 ecological evidence.

    Usage
    -----
    miner = EcologicalMiner(graph)
    miner.run(
        coverage_tsv="coverage_table.tsv",
        annotations_tsv="functional_annotations.tsv",
    )
    # graph now contains co_abundance_guild edges
    """

    def __init__(
        self,
        graph: Optional[nx.MultiGraph] = None,
        r_threshold: float = SPEARMAN_THRESHOLD,
        p_threshold: float = P_VALUE_THRESHOLD,
        pathway_map: Optional[dict[str, str]] = None,
    ):
        self.graph: nx.MultiGraph = graph if graph is not None else nx.MultiGraph()
        self._correlator = CoAbundanceCorrelator(r_threshold, p_threshold)
        self._pathway_checker = PathwayComplementarityChecker(pathway_map)
        self._injector = EcologicalGraphInjector(self.graph)

    def run(
        self,
        coverage_tsv: str | Path,
        annotations_tsv: Optional[str | Path] = None,
    ) -> nx.MultiGraph:
        """
        Run the full Stage 2 ecological evidence pipeline.

        Parameters
        ----------
        coverage_tsv : str | Path
            Coverage table (Contig_ID + sample columns).
        annotations_tsv : str | Path, optional
            Functional annotations table. If omitted, metabolic_match
            defaults to False for all pairs.

        Returns
        -------
        nx.MultiGraph
            Updated graph with co_abundance_guild edges.
        """
        # Module 4 — Correlation
        pairs = self._correlator.correlate(coverage_tsv)

        if not pairs:
            logger.warning("No correlated pairs found above thresholds.")
            return self.graph

        # Module 5 — Pathway complementarity
        if annotations_tsv is not None:
            self._pathway_checker.load_annotations(annotations_tsv)
            pairs = self._pathway_checker.annotate_pairs(pairs)
        else:
            logger.warning(
                "No annotations_tsv provided — metabolic_match set to False for all pairs."
            )
            for p in pairs:
                p["metabolic_match"] = False

        # Module 6 — Graph injection
        self._injector.inject(pairs)

        return self.graph
