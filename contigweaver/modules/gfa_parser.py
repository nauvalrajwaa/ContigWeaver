"""
Module 1: GFA Parser — Physical Evidence Layer
===============================================
Parses a GFA (Graphical Fragment Assembly) file to extract:
  - Segment (S) lines  → contig lengths (node attributes)
  - Link (L) lines     → overlap connections (edges, type='physical_overlap')

Injects results into a NetworkX MultiGraph.
"""

import gzip
import logging
from pathlib import Path
from typing import Optional

import networkx as nx

logger = logging.getLogger(__name__)


class GFAParser:
    """Parse a GFA v1 file and populate a NetworkX graph with physical overlap edges."""

    def __init__(self, graph: Optional[nx.MultiGraph] = None):
        """
        Parameters
        ----------
        graph : nx.MultiGraph, optional
            An existing NetworkX graph to inject data into.
            If None, a new MultiGraph is created.
        """
        self.graph: nx.MultiGraph = graph if graph is not None else nx.MultiGraph()
        self._segment_lengths: dict[str, int] = {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def parse(self, gfa_path: str | Path) -> nx.MultiGraph:
        """
        Parse the GFA file at *gfa_path* and add nodes/edges to ``self.graph``.

        Parameters
        ----------
        gfa_path : str | Path
            Path to the assembly graph GFA file.

        Returns
        -------
        nx.MultiGraph
            The populated graph (same object as ``self.graph``).

        Raises
        ------
        FileNotFoundError
            If *gfa_path* does not exist.
        ValueError
            If the file contains no valid S or L lines.
        """
        gfa_path = Path(gfa_path)
        if not gfa_path.exists():
            raise FileNotFoundError(
                f"GFA file not found: {gfa_path}\n"
                "Please supply the assembly graph produced by MEGAHIT or metaSPAdes."
            )

        logger.info("Parsing GFA file: %s", gfa_path)
        segments_found = 0
        links_found = 0

        open_fn = gzip.open if str(gfa_path).endswith(".gz") else open
        with open_fn(gfa_path, "rt") as fh:
            for lineno, raw_line in enumerate(fh, start=1):
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue

                record_type = line[0]
                fields = line.split("\t")

                if record_type == "S":
                    self._process_segment(fields, lineno)
                    segments_found += 1
                elif record_type == "L":
                    self._process_link(fields, lineno)
                    links_found += 1
                # H (header), P (path), W (walk) lines are intentionally skipped

        if segments_found == 0:
            raise ValueError(
                f"No Segment (S) lines found in {gfa_path}. "
                "Verify the file is a valid GFA v1 assembly graph."
            )

        logger.info(
            "GFA parsing complete — %d segments, %d links ingested.",
            segments_found,
            links_found,
        )
        return self.graph

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _process_segment(self, fields: list[str], lineno: int) -> None:
        """Handle a Segment (S) line.

        GFA v1 S-line format:
            S  <name>  <sequence>  [optional tags …]

        Length is taken from the LN tag if present, otherwise from the
        sequence string itself (using '*' as placeholder for unknown).
        """
        if len(fields) < 3:
            logger.warning("Malformed S line at line %d — skipping.", lineno)
            return

        name = fields[1]
        sequence = fields[2]

        # Try to read the LN:i:<length> optional tag first
        length = self._extract_length_tag(fields[3:])
        if length is None:
            length = len(sequence) if sequence != "*" else 0

        self._segment_lengths[name] = length

        # Add/update node in the graph
        if self.graph.has_node(name):
            self.graph.nodes[name]["length"] = length
        else:
            self.graph.add_node(name, length=length, node_type="unknown")

        logger.debug("Segment %s  length=%d", name, length)

    def _process_link(self, fields: list[str], lineno: int) -> None:
        """Handle a Link (L) line.

        GFA v1 L-line format:
            L  <from>  <from_orient>  <to>  <to_orient>  <overlap_CIGAR>
        """
        if len(fields) < 6:
            logger.warning("Malformed L line at line %d — skipping.", lineno)
            return

        from_node = fields[1]
        from_orient = fields[2]   # '+' or '-'
        to_node = fields[3]
        to_orient = fields[4]     # '+' or '-'
        overlap = fields[5]       # CIGAR string, e.g. '55M'

        # Ensure both endpoint nodes exist (they should already from S lines,
        # but GFA does not guarantee ordering)
        for node in (from_node, to_node):
            if not self.graph.has_node(node):
                length = self._segment_lengths.get(node, 0)
                self.graph.add_node(node, length=length, node_type="unknown")
                logger.debug(
                    "Node %s auto-created from L line (S line may come later).", node
                )

        self.graph.add_edge(
            from_node,
            to_node,
            type="physical_overlap",
            from_orient=from_orient,
            to_orient=to_orient,
            overlap_cigar=overlap,
        )
        logger.debug(
            "Physical overlap: %s(%s) → %s(%s)  CIGAR=%s",
            from_node, from_orient, to_node, to_orient, overlap,
        )

    @staticmethod
    def _extract_length_tag(optional_fields: list[str]) -> Optional[int]:
        """Return the integer value of the LN:i:<n> tag, or None if absent."""
        for field in optional_fields:
            if field.startswith("LN:i:"):
                try:
                    return int(field[5:])
                except ValueError:
                    pass
        return None


# ---------------------------------------------------------------------------
# Convenience function
# ---------------------------------------------------------------------------

def parse_gfa(gfa_path: str | Path, graph: Optional[nx.MultiGraph] = None) -> nx.MultiGraph:
    """
    Parse a GFA file and return (or update) a NetworkX MultiGraph.

    Parameters
    ----------
    gfa_path : str | Path
        Path to the .gfa assembly graph file.
    graph : nx.MultiGraph, optional
        Existing graph to extend; a new one is created if not supplied.

    Returns
    -------
    nx.MultiGraph
        Graph populated with physical overlap edges.
    """
    parser = GFAParser(graph=graph)
    return parser.parse(gfa_path)
