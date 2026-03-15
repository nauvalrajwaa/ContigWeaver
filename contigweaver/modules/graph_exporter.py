"""
Module 3: Graph Integration & Export
=====================================
Merges all evidence layers (physical + CRISPR + ecological) into a single
NetworkX MultiGraph and exports:

  1. ``contigweaver_edges.tsv`` — Tabular edge list for downstream analysis.
  2. ``contigweaver_network.html`` — Interactive PyVis/vis.js visualization.

Visual encoding rules
---------------------
Nodes:
  • Circle  — Bacterial / Unknown contigs
  • Triangle — Viral contigs
  • Size ∝ log(contig_length + 1), scaled to [10, 60] px

Edges:
  • physical_overlap     → solid thick gray (#CCCCCC)
  • crispr_targeting     → dashed red (#FF0000), directed arrow bacteria→virus
  • co_abundance_guild   → solid light blue (#add8e6);
                           thick dark blue (#00008b) when metabolic_match=True
"""

import logging
import math
from pathlib import Path
from typing import Optional

import networkx as nx

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Visual style constants
# ---------------------------------------------------------------------------

EDGE_STYLES: dict[str, dict] = {
    "physical_overlap": {
        "color": "#CCCCCC",
        "width": 4,
        "dashes": False,
        "arrows": "",
        "title_prefix": "Physical Overlap",
    },
    "crispr_targeting": {
        "color": "#FF0000",
        "width": 2,
        "dashes": True,
        "arrows": "to",
        "title_prefix": "CRISPR Targeting",
    },
    "co_abundance_guild": {
        "color": "#add8e6",
        "width": 3,
        "dashes": False,
        "arrows": "",
        "title_prefix": "Co-Abundance",
    },
    "co_abundance_guild_metabolic": {
        "color": "#00008b",
        "width": 5,
        "dashes": False,
        "arrows": "",
        "title_prefix": "Co-Abundance (metabolic match)",
    },
}

NODE_SIZE_MIN = 10
NODE_SIZE_MAX = 60


# ===========================================================================
# Graph Exporter
# ===========================================================================

class GraphExporter:
    """Export a NetworkX MultiGraph to TSV and interactive HTML."""

    def __init__(self, graph: nx.MultiGraph):
        self.graph = graph

    # ------------------------------------------------------------------
    # TSV Export
    # ------------------------------------------------------------------

    def export_tsv(self, output_path: str | Path) -> Path:
        """
        Write the edge list to a TSV file.

        Columns:
            Source_Contig, Target_Contig, Evidence_Type,
            Weight_or_Identity, Source_Length, Target_Length

        Returns
        -------
        Path
            Path to the written TSV file.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        rows: list[str] = [
            "Source_Contig\tTarget_Contig\tEvidence_Type\t"
            "Weight_or_Identity\tSource_Length\tTarget_Length"
        ]

        for u, v, data in self.graph.edges(data=True):
            edge_type = data.get("type", "unknown")
            evidence_type = self._map_evidence_type(edge_type, data)
            weight = self._get_weight(edge_type, data)
            src_len = self.graph.nodes[u].get("length", 0)
            tgt_len = self.graph.nodes[v].get("length", 0)
            rows.append(
                f"{u}\t{v}\t{evidence_type}\t{weight}\t{src_len}\t{tgt_len}"
            )

        output_path.write_text("\n".join(rows) + "\n")
        logger.info("TSV edge list written to %s (%d edges).", output_path, len(rows) - 1)
        return output_path

    @staticmethod
    def _map_evidence_type(edge_type: str, data: dict) -> str:
        if edge_type == "physical_overlap":
            return "Physical"
        if edge_type == "crispr_targeting":
            return "CRISPR"
        if edge_type == "co_abundance_guild":
            return "Co-Abundance"
        return edge_type

    @staticmethod
    def _get_weight(edge_type: str, data: dict) -> str:
        if edge_type == "physical_overlap":
            return data.get("overlap_cigar", "NA")
        if edge_type == "crispr_targeting":
            return f"{data.get('identity', 0.0):.2f}"
        if edge_type == "co_abundance_guild":
            return f"{data.get('weight', 0.0):.4f}"
        return "NA"

    # ------------------------------------------------------------------
    # HTML / PyVis Export
    # ------------------------------------------------------------------

    def export_html(
        self,
        output_path: str | Path,
        viral_contigs: Optional[set[str]] = None,
        height: str = "800px",
        width: str = "100%",
    ) -> Path:
        """
        Render the graph as an interactive HTML file using PyVis.

        Parameters
        ----------
        output_path : str | Path
            Where to save the HTML file.
        viral_contigs : set[str], optional
            Set of contig IDs known to be viral (used for triangle shape).
            Falls back to ``node_type=='viral'`` attribute if not supplied.
        height : str
            Canvas height (CSS value).
        width : str
            Canvas width (CSS value).

        Returns
        -------
        Path
            Path to the written HTML file.

        Raises
        ------
        ImportError
            If pyvis is not installed.
        """
        try:
            from pyvis.network import Network  # type: ignore
        except ImportError as exc:
            raise ImportError(
                "PyVis is required for HTML export.\n"
                "Install it with: pip install pyvis"
            ) from exc

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        viral_set = viral_contigs or set()

        net = Network(
            height=height,
            width=width,
            directed=True,
            bgcolor="#1a1a2e",
            font_color="white",
            notebook=False,
        )
        net.force_atlas_2based()

        # Add nodes
        for node, attrs in self.graph.nodes(data=True):
            node_type = attrs.get("node_type", "unknown")
            is_viral = (node in viral_set) or (node_type == "viral")
            length = attrs.get("length", 1)
            size = self._length_to_size(length)

            shape = "triangle" if is_viral else "dot"
            color = "#e74c3c" if is_viral else "#3498db"
            title = (
                f"<b>{node}</b><br>"
                f"Type: {'Viral' if is_viral else 'Bacterial/Unknown'}<br>"
                f"Length: {length:,} bp"
            )

            net.add_node(
                node,
                label=node,
                shape=shape,
                color=color,
                size=size,
                title=title,
            )

        # Add edges
        for u, v, data in self.graph.edges(data=True):
            edge_type = data.get("type", "unknown")
            style = self._get_edge_style(edge_type, data)
            title = self._build_edge_tooltip(edge_type, data)

            net.add_edge(
                u, v,
                color=style["color"],
                width=style["width"],
                dashes=style["dashes"],
                arrows=style["arrows"],
                title=title,
            )

        # Configure physics & interaction
        net.set_options("""
        {
          "physics": {
            "forceAtlas2Based": {
              "gravitationalConstant": -50,
              "centralGravity": 0.01,
              "springLength": 120,
              "springConstant": 0.08
            },
            "solver": "forceAtlas2Based",
            "stabilization": { "iterations": 150 }
          },
          "interaction": {
            "hover": true,
            "tooltipDelay": 100,
            "navigationButtons": true,
            "keyboard": true
          },
          "edges": {
            "smooth": { "type": "continuous" }
          }
        }
        """)

        net.write_html(str(output_path))
        logger.info("Interactive HTML visualization written to %s.", output_path)
        return output_path

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _length_to_size(length: int) -> float:
        """Map contig length to node size on a log scale."""
        if length <= 0:
            return NODE_SIZE_MIN
        log_val = math.log1p(length)
        # Approximate typical range: 1 bp → size 10, 1 Mbp → size 60
        log_min, log_max = math.log1p(1), math.log1p(1_000_000)
        scaled = (log_val - log_min) / (log_max - log_min)
        return NODE_SIZE_MIN + scaled * (NODE_SIZE_MAX - NODE_SIZE_MIN)

    @staticmethod
    def _get_edge_style(edge_type: str, data: dict) -> dict:
        if edge_type == "co_abundance_guild" and data.get("metabolic_match"):
            return EDGE_STYLES["co_abundance_guild_metabolic"]
        return EDGE_STYLES.get(edge_type, {
            "color": "#888888",
            "width": 1,
            "dashes": False,
            "arrows": "",
            "title_prefix": edge_type,
        })

    @staticmethod
    def _build_edge_tooltip(edge_type: str, data: dict) -> str:
        if edge_type == "physical_overlap":
            return (
                f"<b>Physical Overlap</b><br>"
                f"CIGAR: {data.get('overlap_cigar', 'NA')}<br>"
                f"Orientation: {data.get('from_orient', '?')}/{data.get('to_orient', '?')}"
            )
        if edge_type == "crispr_targeting":
            return (
                f"<b>CRISPR Targeting</b><br>"
                f"Identity: {data.get('identity', 0.0):.2f}%<br>"
                f"Coverage: {data.get('coverage', 0.0):.2%}<br>"
                f"E-value: {data.get('evalue', 'NA')}"
            )
        if edge_type == "co_abundance_guild":
            tooltip = (
                f"<b>Co-Abundance Guild</b><br>"
                f"Spearman r: {data.get('weight', 0.0):.4f}<br>"
                f"Metabolic match: {data.get('metabolic_match', False)}"
            )
            return tooltip
        return f"<b>{edge_type}</b>"


# ---------------------------------------------------------------------------
# Convenience function
# ---------------------------------------------------------------------------

def export_graph(
    graph: nx.MultiGraph,
    tsv_path: str | Path = "contigweaver_edges.tsv",
    html_path: str | Path = "contigweaver_network.html",
    viral_contigs: Optional[set[str]] = None,
) -> tuple[Path, Path]:
    """
    Export a ContigNexus graph to TSV and interactive HTML.

    Parameters
    ----------
    graph : nx.MultiGraph
        The fully populated interaction graph.
    tsv_path : str | Path
        Output path for the TSV edge list.
    html_path : str | Path
        Output path for the HTML visualization.
    viral_contigs : set[str], optional
        Explicit set of viral contig IDs for node shaping.

    Returns
    -------
    tuple[Path, Path]
        (tsv_path, html_path)
    """
    exporter = GraphExporter(graph)
    tsv = exporter.export_tsv(tsv_path)
    html = exporter.export_html(html_path, viral_contigs=viral_contigs)
    return tsv, html
