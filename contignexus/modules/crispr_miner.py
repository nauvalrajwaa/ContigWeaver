"""
Module 2: CRISPR-Phage Miner — Biological Evidence Layer
=========================================================
Workflow:
  1. Run ``minced`` on ``contigs.fasta`` to detect CRISPR arrays.
  2. Parse minced output → extract spacer sequences → write ``spacers.fasta``.
  3. Run ``blastn`` (blastn-short) to align spacers against ``viral_contigs.fasta``.
  4. Filter BLAST hits: identity > 95 %, near-full coverage (>= 90 % of spacer length).
  5. Inject filtered pairs into a NetworkX graph as edges with
     type='crispr_targeting', from bacterial contig → viral contig.
"""

import logging
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import networkx as nx

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BLAST_IDENTITY_THRESHOLD = 95.0   # percent identity
BLAST_COVERAGE_THRESHOLD = 0.90   # fraction of query (spacer) length aligned


# ===========================================================================
# MinCED Wrapper
# ===========================================================================

class MincedRunner:
    """Lightweight wrapper around the ``minced`` CRISPR detection tool."""

    def __init__(self, minced_bin: str = "minced"):
        """
        Parameters
        ----------
        minced_bin : str
            Path or name of the minced executable.
        """
        self.minced_bin = minced_bin

    def run(self, contigs_fasta: str | Path, output_dir: str | Path) -> Path:
        """
        Run minced on *contigs_fasta* and write output to *output_dir*.

        Returns
        -------
        Path
            Path to the minced ``.txt`` output file.

        Raises
        ------
        FileNotFoundError
            If *contigs_fasta* is missing.
        RuntimeError
            If minced exits with a non-zero return code.
        """
        contigs_fasta = Path(contigs_fasta)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if not contigs_fasta.exists():
            raise FileNotFoundError(
                f"Contigs FASTA not found: {contigs_fasta}"
            )

        out_txt = output_dir / "minced_output.txt"
        out_gff = output_dir / "minced_output.gff"

        cmd = [
            self.minced_bin,
            str(contigs_fasta),
            str(out_txt),
            str(out_gff),
        ]
        logger.info("Running minced: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(
                f"minced failed (exit {result.returncode}):\n"
                f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

        logger.info("minced finished. Output: %s", out_txt)
        return out_txt


# ===========================================================================
# Spacer Extractor
# ===========================================================================

class SpacerExtractor:
    """Parse minced output text and extract spacer sequences."""

    # Example minced line with spacers:
    #   Spacers: ACGT... TTGA... ...
    # Minced text output format repeats CRISPR array blocks per contig.
    _SPACER_RE = re.compile(r"^\s+\d+\s+([ACGTacgtNn]+)\s+\[", re.MULTILINE)
    _CONTIG_RE = re.compile(r"^Sequence '(.+?)' \(", re.MULTILINE)

    def extract_to_fasta(
        self,
        minced_txt: str | Path,
        spacers_fasta: str | Path,
    ) -> int:
        """
        Parse *minced_txt* and write all spacers to *spacers_fasta*.

        Returns
        -------
        int
            Number of spacers extracted.

        Raises
        ------
        FileNotFoundError
            If *minced_txt* does not exist.
        """
        minced_txt = Path(minced_txt)
        spacers_fasta = Path(spacers_fasta)

        if not minced_txt.exists():
            raise FileNotFoundError(f"minced output not found: {minced_txt}")

        text = minced_txt.read_text()
        spacer_count = 0
        current_contig = "unknown"

        with spacers_fasta.open("w") as out:
            for line in text.splitlines():
                # Detect contig header lines
                contig_match = self._CONTIG_RE.match(line)
                if contig_match:
                    current_contig = contig_match.group(1)
                    continue

                # Detect spacer sequence lines:
                # Minced formats each repeat-spacer row as:
                #   <idx>  <repeat>  [ <length>, <length> ]  <spacer>
                # We split and take the last token before the bracket as spacer
                spacer_match = self._parse_spacer_line(line)
                if spacer_match:
                    spacer_count += 1
                    header = f">spacer_{spacer_count}|contig={current_contig}"
                    out.write(f"{header}\n{spacer_match}\n")

        logger.info("Extracted %d spacers to %s", spacer_count, spacers_fasta)
        return spacer_count

    @staticmethod
    def _parse_spacer_line(line: str) -> Optional[str]:
        """
        Extract spacer sequence from a minced body line.

        Minced output rows look like::
            1   GCCTCCCACTTATCCGGATGATCAT   [ 25, 36 ]   AGTTGGAAGATGCTTTTAAGAAATCAGAAT

        The *spacer* is the second DNA token (after the repeat).
        Returns None if the line does not match.
        """
        stripped = line.strip()
        if not stripped:
            return None

        # Must start with a digit (row index)
        parts = stripped.split()
        if not parts or not parts[0].isdigit():
            return None

        # Collect DNA tokens (only [ACGTacgtNn]+)
        dna_re = re.compile(r"^[ACGTacgtNn]+$")
        dna_tokens = [p for p in parts if dna_re.match(p)]

        # Row: index  repeat  [len, len]  spacer
        # dna_tokens[0] = repeat, dna_tokens[1] = spacer (if present)
        if len(dna_tokens) >= 2:
            return dna_tokens[1].upper()
        return None


# ===========================================================================
# BLAST Runner
# ===========================================================================

class BlastRunner:
    """Wrapper to run blastn (blastn-short) and parse results."""

    def __init__(
        self,
        blastn_bin: str = "blastn",
        makeblastdb_bin: str = "makeblastdb",
    ):
        self.blastn_bin = blastn_bin
        self.makeblastdb_bin = makeblastdb_bin

    def build_db(self, viral_fasta: str | Path, db_dir: str | Path) -> Path:
        """
        Build a BLAST nucleotide database from *viral_fasta*.

        Returns
        -------
        Path
            Path to the database prefix (suitable for -db argument).
        """
        viral_fasta = Path(viral_fasta)
        db_dir = Path(db_dir)
        db_dir.mkdir(parents=True, exist_ok=True)

        if not viral_fasta.exists():
            raise FileNotFoundError(
                f"Viral contigs FASTA not found: {viral_fasta}"
            )

        db_path = db_dir / "viral_db"
        cmd = [
            self.makeblastdb_bin,
            "-in", str(viral_fasta),
            "-dbtype", "nucl",
            "-out", str(db_path),
            "-title", "viral_contigs",
        ]
        logger.info("Building BLAST DB: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"makeblastdb failed (exit {result.returncode}):\n{result.stderr}"
            )
        return db_path

    def run_blastn_short(
        self,
        query_fasta: str | Path,
        db_path: str | Path,
        out_tsv: str | Path,
        evalue: float = 1e-3,
    ) -> Path:
        """
        Run blastn-short on *query_fasta* against *db_path*.

        Output format: qseqid sseqid pident length qlen slen evalue bitscore

        Returns
        -------
        Path
            Path to the tabular BLAST output file.
        """
        query_fasta = Path(query_fasta)
        out_tsv = Path(out_tsv)

        cmd = [
            self.blastn_bin,
            "-task", "blastn-short",
            "-query", str(query_fasta),
            "-db", str(db_path),
            "-out", str(out_tsv),
            "-evalue", str(evalue),
            "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
            "-perc_identity", str(BLAST_IDENTITY_THRESHOLD),
        ]
        logger.info("Running BLAST: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"blastn failed (exit {result.returncode}):\n{result.stderr}"
            )
        logger.info("BLAST finished. Hits file: %s", out_tsv)
        return out_tsv


# ===========================================================================
# BLAST Result Filter
# ===========================================================================

class BlastFilter:
    """Filter BLAST tabular results by identity and coverage thresholds."""

    def filter(
        self,
        blast_tsv: str | Path,
        identity_threshold: float = BLAST_IDENTITY_THRESHOLD,
        coverage_threshold: float = BLAST_COVERAGE_THRESHOLD,
    ) -> list[dict]:
        """
        Read *blast_tsv* and return rows passing both filters.

        Parameters
        ----------
        blast_tsv : str | Path
            Tab-separated BLAST output (outfmt 6 with columns:
            qseqid sseqid pident length qlen slen evalue bitscore).
        identity_threshold : float
            Minimum percent identity (default 95.0).
        coverage_threshold : float
            Minimum fraction of query length aligned (default 0.90).

        Returns
        -------
        list[dict]
            Filtered records as dictionaries.
        """
        blast_tsv = Path(blast_tsv)
        if not blast_tsv.exists():
            logger.warning("BLAST output not found: %s — returning empty list.", blast_tsv)
            return []

        results = []
        with blast_tsv.open() as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 8:
                    continue

                try:
                    pident = float(parts[2])
                    aln_len = int(parts[3])
                    qlen = int(parts[4])
                except ValueError:
                    continue

                coverage = aln_len / qlen if qlen > 0 else 0.0

                if pident >= identity_threshold and coverage >= coverage_threshold:
                    results.append(
                        {
                            "qseqid": parts[0],
                            "sseqid": parts[1],
                            "pident": pident,
                            "aln_len": aln_len,
                            "qlen": qlen,
                            "slen": int(parts[5]),
                            "evalue": float(parts[6]),
                            "bitscore": float(parts[7]),
                            "coverage": coverage,
                        }
                    )

        logger.info(
            "BLAST filter: %d hits passed (identity>=%.1f%%, coverage>=%.0f%%).",
            len(results), identity_threshold, coverage_threshold * 100,
        )
        return results


# ===========================================================================
# High-Level CRISPR-Phage Miner
# ===========================================================================

class CRISPRPhageMiner:
    """
    Orchestrates the full CRISPR-Phage evidence pipeline:
      minced → spacer extraction → BLAST → filter → NetworkX injection.
    """

    def __init__(
        self,
        graph: Optional[nx.MultiGraph] = None,
        minced_bin: str = "minced",
        blastn_bin: str = "blastn",
        makeblastdb_bin: str = "makeblastdb",
    ):
        self.graph: nx.MultiGraph = graph if graph is not None else nx.MultiGraph()
        self._minced = MincedRunner(minced_bin)
        self._spacer_extractor = SpacerExtractor()
        self._blast = BlastRunner(blastn_bin, makeblastdb_bin)
        self._filter = BlastFilter()

    def run(
        self,
        contigs_fasta: str | Path,
        viral_contigs_fasta: str | Path,
        work_dir: str | Path = "contignexus_workdir",
    ) -> nx.MultiGraph:
        """
        Execute the complete CRISPR-Phage mining pipeline.

        Parameters
        ----------
        contigs_fasta : str | Path
            All assembled contigs (input to minced).
        viral_contigs_fasta : str | Path
            Viral contig subset (BLAST database subject).
        work_dir : str | Path
            Working directory for intermediate files.

        Returns
        -------
        nx.MultiGraph
            Updated graph with crispr_targeting edges added.
        """
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        # Step 1 — Run minced
        minced_out = self._minced.run(contigs_fasta, work_dir / "minced")

        # Step 2 — Extract spacers to FASTA
        spacers_fasta = work_dir / "spacers.fasta"
        n_spacers = self._spacer_extractor.extract_to_fasta(minced_out, spacers_fasta)

        if n_spacers == 0:
            logger.warning(
                "No CRISPR spacers found in %s. "
                "Skipping BLAST step — no crispr_targeting edges will be added.",
                contigs_fasta,
            )
            return self.graph

        # Step 3 — Build BLAST DB & run blastn-short
        db_path = self._blast.build_db(viral_contigs_fasta, work_dir / "blastdb")
        blast_out = work_dir / "spacer_vs_viral.tsv"
        self._blast.run_blastn_short(spacers_fasta, db_path, blast_out)

        # Step 4 — Filter results
        hits = self._filter.filter(blast_out)

        # Step 5 — Inject edges into graph
        self._inject_edges(hits)

        return self.graph

    def _inject_edges(self, hits: list[dict]) -> None:
        """
        Add crispr_targeting edges to the graph.

        The edge source is the *bacterial* contig (spacer origin), and the
        target is the *viral* contig (BLAST subject). The spacer header
        carries ``contig=<name>`` which is used to resolve the bacterial side.
        """
        added = 0
        for hit in hits:
            # spacer header format: spacer_<n>|contig=<bacterial_contig>
            bact_contig = self._parse_contig_from_spacer_id(hit["qseqid"])
            viral_contig = hit["sseqid"]

            if not self.graph.has_node(bact_contig):
                self.graph.add_node(bact_contig, length=0, node_type="unknown")
            if not self.graph.has_node(viral_contig):
                self.graph.add_node(viral_contig, length=0, node_type="viral")
            else:
                # Mark target as viral
                self.graph.nodes[viral_contig]["node_type"] = "viral"

            self.graph.add_edge(
                bact_contig,
                viral_contig,
                type="crispr_targeting",
                identity=hit["pident"],
                coverage=hit["coverage"],
                evalue=hit["evalue"],
            )
            added += 1

        logger.info("Injected %d crispr_targeting edges into graph.", added)

    @staticmethod
    def _parse_contig_from_spacer_id(spacer_id: str) -> str:
        """Extract the bacterial contig name from a spacer FASTA header."""
        # Header format: spacer_<n>|contig=<name>
        if "|contig=" in spacer_id:
            return spacer_id.split("|contig=")[-1]
        return spacer_id
