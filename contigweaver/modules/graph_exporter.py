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
    "segment_membership": {
        "color": "#f39c12",
        "width": 2,
        "dashes": True,
        "arrows": "",
        "title_prefix": "Segment Membership",
    },
}

NODE_SIZE_MIN = 10
NODE_SIZE_MAX = 60
HTML_MAX_NODES = 1500
HTML_MAX_EDGES = 4000


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
        if edge_type == "segment_membership":
            return "Segment-Membership"
        return edge_type

    @staticmethod
    def _get_weight(edge_type: str, data: dict) -> str:
        if edge_type == "physical_overlap":
            return data.get("overlap_cigar", "NA")
        if edge_type == "crispr_targeting":
            return f"{data.get('identity', 0.0):.2f}"
        if edge_type == "co_abundance_guild":
            return f"{data.get('weight', 0.0):.4f}"
        if edge_type == "segment_membership":
            return str(data.get("anchor_length", "NA"))
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
        html_graph, html_meta = self._prepare_html_graph()

        net = Network(
            height=height,
            width=width,
            directed=True,
            bgcolor="#1a1a2e",
            notebook=False,
        )
        net.force_atlas_2based()
        if html_meta["sampled"]:
            net.heading = (
                "ContigWeaver Focused HTML View "
                f"({html_meta['nodes']} nodes, {html_meta['edges']} edges from "
                f"{self.graph.number_of_nodes():,}/{self.graph.number_of_edges():,})"
            )

        # Add nodes
        for node, attrs in html_graph.nodes(data=True):
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
        for u, v, data in html_graph.edges(data=True):
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

    def _prepare_html_graph(
        self,
        max_nodes: int = HTML_MAX_NODES,
        max_edges: int = HTML_MAX_EDGES,
        focus_hops: int = 1,
    ) -> tuple[nx.MultiGraph, dict[str, int | bool]]:
        if (
            self.graph.number_of_nodes() <= max_nodes
            and self.graph.number_of_edges() <= max_edges
        ):
            return self.graph, {
                "sampled": False,
                "nodes": self.graph.number_of_nodes(),
                "edges": self.graph.number_of_edges(),
            }

        focus_nodes = self._focus_nodes()
        if not focus_nodes:
            focus_nodes = self._top_degree_nodes(limit=max_nodes)

        selected_nodes: set[str] = set(focus_nodes)
        frontier = list(focus_nodes)
        for _ in range(max(focus_hops, 0)):
            next_frontier: list[str] = []
            for node in frontier:
                for raw_neighbor in self.graph.neighbors(node):
                    neighbor = str(raw_neighbor)
                    if len(selected_nodes) >= max_nodes:
                        break
                    if neighbor in selected_nodes:
                        continue
                    selected_nodes.add(neighbor)
                    next_frontier.append(neighbor)
                if len(selected_nodes) >= max_nodes:
                    break
            frontier = next_frontier
            if not frontier or len(selected_nodes) >= max_nodes:
                break

        if len(selected_nodes) < max_nodes:
            for node in self._top_degree_nodes(limit=max_nodes * 2):
                if len(selected_nodes) >= max_nodes:
                    break
                selected_nodes.add(node)

        subgraph = nx.MultiGraph()
        for node in selected_nodes:
            if self.graph.has_node(node):
                subgraph.add_node(node, **self.graph.nodes[node])
        for u, v, key, data in self.graph.edges(keys=True, data=True):
            if u in selected_nodes and v in selected_nodes:
                subgraph.add_edge(u, v, key=key, **data)
        if subgraph.number_of_edges() <= max_edges:
            return subgraph, {
                "sampled": True,
                "nodes": subgraph.number_of_nodes(),
                "edges": subgraph.number_of_edges(),
            }

        edge_rows: list[tuple[tuple[str, str, int], dict, tuple[int, int, int]]] = []
        degrees = {str(node): int(degree) for node, degree in self.graph.degree()}
        for u, v, key, data in subgraph.edges(keys=True, data=True):
            edge_type = data.get("type", "unknown")
            priority = 0 if edge_type != "physical_overlap" else 1
            bridge_priority = 0 if edge_type == "segment_membership" else 1
            score = -(degrees.get(u, 0) + degrees.get(v, 0))
            edge_rows.append(((u, v, key), data, (priority, bridge_priority, score)))

        edge_rows.sort(key=lambda item: item[2])
        limited = nx.MultiGraph()
        for node, attrs in subgraph.nodes(data=True):
            limited.add_node(node, **attrs)
        for (u, v, _), data, _ in edge_rows[:max_edges]:
            limited.add_edge(u, v, **data)

        return limited, {
            "sampled": True,
            "nodes": limited.number_of_nodes(),
            "edges": limited.number_of_edges(),
        }

    def _focus_nodes(self) -> list[str]:
        focus_nodes: set[str] = set()
        for u, v, data in self.graph.edges(data=True):
            edge_type = data.get("type")
            if edge_type and edge_type != "physical_overlap":
                focus_nodes.add(u)
                focus_nodes.add(v)
        return sorted(focus_nodes)

    def _top_degree_nodes(self, limit: int) -> list[str]:
        ranked = sorted(
            ((str(node), int(degree)) for node, degree in self.graph.degree()),
            key=lambda item: item[1],
            reverse=True,
        )
        return [node for node, _ in ranked[:limit]]

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
        if edge_type == "segment_membership":
            return (
                f"<b>Segment Membership</b><br>"
                f"Anchor length: {data.get('anchor_length', 'NA')} bp<br>"
                f"Match type: {data.get('match_type', 'NA')}<br>"
                f"Orientation: {data.get('orientation', '?')}"
            )
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
