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
import hashlib
from collections import Counter
from datetime import datetime
from html import escape
from os.path import relpath
from pathlib import Path
from typing import Mapping, Optional, Sequence

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
    "taxonomic_match": {
        "color": "#2ecc71",
        "width": 2,
        "dashes": True,
        "arrows": "",
        "title_prefix": "Taxonomic Match",
    },
    "functional_operon": {
        "color": "#ff7f0e",
        "width": 3,
        "dashes": False,
        "arrows": "to",
        "title_prefix": "Functional Operon",
    },
    "bin_membership": {
        "color": "#7f8c8d",
        "width": 1,
        "dashes": True,
        "arrows": "",
        "title_prefix": "Bin Membership",
    },
}

NODE_SIZE_MIN = 10
NODE_SIZE_MAX = 60
HTML_MAX_NODES = 1500
HTML_MAX_EDGES = 4000
BIN_COLOR_PALETTE = (
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
)


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
        if edge_type == "taxonomic_match":
            return "Taxonomic-Match"
        if edge_type == "functional_operon":
            return "Functional-Operon"
        if edge_type == "bin_membership":
            return "Bin-Membership"
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
        if edge_type == "taxonomic_match":
            return f"{data.get('confidence_min', 0.0):.3f}"
        if edge_type == "functional_operon":
            return f"{data.get('score', 0.0):.3f}"
        if edge_type == "bin_membership":
            return str(data.get("bin_id", "NA"))
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
            color = self._node_color(str(node), attrs, is_viral)
            title = self._build_node_tooltip(str(node), attrs, is_viral)

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
                edge_type=edge_type,
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
        self._inject_edge_type_filters(
            output_path,
            sorted({str(data.get("type", "unknown")) for _, _, data in html_graph.edges(data=True)}),
        )
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
            if edge_type in {
                "crispr_targeting",
                "co_abundance_guild",
                "functional_operon",
                "taxonomic_match",
                "segment_membership",
            }:
                priority = 0
            elif edge_type == "bin_membership":
                priority = 1
            else:
                priority = 2

            if edge_type == "segment_membership":
                bridge_priority = 0
            elif edge_type in {
                "crispr_targeting",
                "co_abundance_guild",
                "functional_operon",
                "taxonomic_match",
            }:
                bridge_priority = 1
            elif edge_type == "bin_membership":
                bridge_priority = 2
            else:
                bridge_priority = 3
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
            if edge_type and edge_type not in {"physical_overlap", "bin_membership"}:
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
        if edge_type == "taxonomic_match":
            return (
                f"<b>Taxonomic Match</b><br>"
                f"Rank: {data.get('taxonomy_rank', 'NA')}<br>"
                f"Label: {data.get('taxonomy_label', 'NA')}<br>"
                f"Min confidence: {data.get('confidence_min', 0.0):.3f}"
            )
        if edge_type == "functional_operon":
            return (
                f"<b>Functional Operon</b><br>"
                f"Cas genes: {data.get('cas_genes', 'NA')}<br>"
                f"Support: {data.get('support_mode', 'NA')}<br>"
                f"Distance (hops): {data.get('distance_hops', 'NA')}"
            )
        if edge_type == "bin_membership":
            return (
                f"<b>Bin Membership</b><br>"
                f"Bin: {data.get('bin_id', 'NA')}<br>"
                f"Assignment source: {data.get('assignment_source', 'NA')}<br>"
                f"Rescue support: {data.get('rescue_support_count', 'NA')}"
            )
        return f"<b>{edge_type}</b>"

    @staticmethod
    def _node_color(node_id: str, attrs: dict, is_viral: bool) -> str:
        if is_viral:
            return "#e74c3c"
        bin_id = attrs.get("bin_id")
        if not bin_id:
            return "#3498db"
        digest = hashlib.md5(str(bin_id).encode("utf-8"), usedforsecurity=False).hexdigest()
        index = int(digest[:8], 16) % len(BIN_COLOR_PALETTE)
        return BIN_COLOR_PALETTE[index]

    @staticmethod
    def _build_node_tooltip(node_id: str, attrs: dict, is_viral: bool) -> str:
        length = int(attrs.get("length", 0) or 0)
        lines = [
            f"<b>{escape(node_id)}</b>",
            f"Type: {'Viral' if is_viral else 'Bacterial/Unknown'}",
            f"Length: {length:,} bp",
        ]

        taxonomy_label = attrs.get("taxonomy_label")
        taxonomy_rank = attrs.get("taxonomy_rank")
        if taxonomy_label:
            lines.append(
                "Taxonomy: "
                + escape(str(taxonomy_label))
                + (f" ({escape(str(taxonomy_rank))})" if taxonomy_rank else "")
            )

        confidence = attrs.get("taxonomy_confidence")
        if confidence not in (None, ""):
            try:
                lines.append(f"Taxonomy confidence: {float(confidence):.3f}")
            except (TypeError, ValueError):
                lines.append(f"Taxonomy confidence: {escape(str(confidence))}")

        functional_terms = attrs.get("functional_terms")
        if functional_terms:
            if isinstance(functional_terms, list):
                rendered = ", ".join(str(item) for item in functional_terms[:8])
            else:
                rendered = str(functional_terms)
            lines.append(f"Functions: {escape(rendered)}")

        cas_genes = attrs.get("cas_genes")
        if cas_genes:
            if isinstance(cas_genes, list):
                rendered = ", ".join(str(item) for item in cas_genes[:8])
            else:
                rendered = str(cas_genes)
            lines.append(f"Cas genes: {escape(rendered)}")

        bin_id = attrs.get("bin_id")
        if bin_id:
            lines.append(f"Bin: {escape(str(bin_id))}")
            lines.append(f"Bin status: {escape(str(attrs.get('bin_status', 'binned')))}")
            if attrs.get("rescued"):
                lines.append("Rescued: yes")

        if attrs.get("has_crispr"):
            spacer_count = attrs.get("spacer_count", 0)
            lines.append(f"CRISPR support: yes (spacers={escape(str(spacer_count))})")

        return "<br>".join(lines)

    @staticmethod
    def _edge_type_display_name(edge_type: str) -> str:
        labels = {
            "physical_overlap": "Physical Overlap",
            "crispr_targeting": "CRISPR Targeting",
            "co_abundance_guild": "Co-Abundance",
            "segment_membership": "Segment Membership",
            "taxonomic_match": "Taxonomic Match",
            "functional_operon": "Functional Operon",
            "bin_membership": "Bin Membership",
        }
        return labels.get(edge_type, edge_type)

    @staticmethod
    def _inject_edge_type_filters(output_path: Path, edge_types: list[str]) -> None:
        if not edge_types:
            return

        html = output_path.read_text()
        control_items = "".join(
            (
                '<label style="display:inline-flex;align-items:center;gap:6px;margin-right:12px;">'
                f'<input type="checkbox" data-edge-type="{escape(edge_type)}" checked>'
                f"{escape(GraphExporter._edge_type_display_name(edge_type))}"
                "</label>"
            )
            for edge_type in edge_types
        )
        controls = (
            '<div class="cw-edge-filter" '
            'style="position:fixed;top:10px;left:10px;z-index:9999;background:rgba(255,255,255,0.96);padding:10px 12px;border:1px solid #ccd2d7;border-radius:10px;font-family:Arial,sans-serif;font-size:12px;max-width:70vw;">'
            '<div style="font-weight:600;margin-bottom:6px;">Edge Filters</div>'
            f"{control_items}"
            "</div>"
        )
        script = (
            "<script>(function(){"
            "if(typeof edges==='undefined'){return;}"
            "const controls=Array.from(document.querySelectorAll('.cw-edge-filter input[data-edge-type]'));"
            "if(!controls.length){return;}"
            "const initial={};"
            "edges.get().forEach(function(edge){"
            "initial[edge.id]={edge_type:String(edge.edge_type||'unknown')};"
            "});"
            "function applyFilters(){"
            "const enabled=new Set(controls.filter(function(input){return input.checked;}).map(function(input){return input.getAttribute('data-edge-type')||'';}));"
            "const updates=[];"
            "edges.get().forEach(function(edge){"
            "const meta=initial[edge.id];"
            "const edgeType=(meta&&meta.edge_type)?meta.edge_type:'unknown';"
            "updates.push({id:edge.id,hidden:!enabled.has(edgeType)});"
            "});"
            "edges.update(updates);"
            "if(typeof network!=='undefined'){network.fit({animation:false});}"
            "}"
            "controls.forEach(function(input){input.addEventListener('change',applyFilters);});"
            "applyFilters();"
            "})();</script>"
        )

        if "<body>" in html:
            html = html.replace("<body>", "<body>\n" + controls, 1)
        else:
            html = controls + html

        if "</body>" in html:
            html = html.replace("</body>", script + "\n</body>", 1)
        else:
            html += script

        output_path.write_text(html)


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


def export_index_report(
    graph: nx.MultiGraph,
    output_path: str | Path,
    network_html_path: str | Path,
    edges_tsv_path: str | Path,
    viral_contigs: Optional[set[str]] = None,
    input_paths: Optional[Mapping[str, str | Path | None]] = None,
    related_html_paths: Optional[Sequence[str | Path]] = None,
) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    network_html_path = Path(network_html_path)
    edges_tsv_path = Path(edges_tsv_path)
    related_html_paths = [Path(path) for path in (related_html_paths or [])]

    viral_set = viral_contigs or set()
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    node_count = graph.number_of_nodes()
    edge_count = graph.number_of_edges()
    viral_count = sum(
        1
        for node, attrs in graph.nodes(data=True)
        if node in viral_set or attrs.get("node_type") == "viral"
    )
    evidence_counts = Counter(
        data.get("type", "unknown") for _, _, data in graph.edges(data=True)
    )
    top_nodes = sorted(
        ((str(node), int(degree)) for node, degree in graph.degree()),
        key=lambda item: item[1],
        reverse=True,
    )[:10]

    rows = []
    for source, target, data in list(graph.edges(data=True))[:15]:
        rows.append(
            "<tr>"
            f"<td>{escape(str(source))}</td>"
            f"<td>{escape(str(target))}</td>"
            f"<td>{escape(str(data.get('type', 'unknown')))}</td>"
            f"<td>{escape(str(GraphExporter._get_weight(data.get('type', 'unknown'), data)))}</td>"
            "</tr>"
        )
    sample_rows = "\n".join(rows) or (
        '<tr><td colspan="4">No edges available.</td></tr>'
    )

    input_items = []
    for label, raw_path in (input_paths or {}).items():
        if raw_path is None:
            continue
        path = Path(raw_path)
        path_label = _relative_link(path, output_path.parent)
        status = "present" if path.exists() else "missing"
        input_items.append(
            "<tr>"
            f"<td>{escape(label)}</td>"
            f"<td><a href=\"{escape(path_label)}\">{escape(str(path))}</a></td>"
            f"<td>{status}</td>"
            "</tr>"
        )
    input_table = "\n".join(input_items) or (
        '<tr><td colspan="3">No run inputs recorded.</td></tr>'
    )

    evidence_rows = "\n".join(
        "<tr>"
        f"<td>{escape(edge_type)}</td>"
        f"<td>{count:,}</td>"
        "</tr>"
        for edge_type, count in sorted(evidence_counts.items())
    ) or '<tr><td colspan="2">No evidence edges exported.</td></tr>'

    degree_rows = "\n".join(
        "<tr>"
        f"<td>{escape(node)}</td>"
        f"<td>{degree:,}</td>"
        "</tr>"
        for node, degree in top_nodes
    ) or '<tr><td colspan="2">No nodes available.</td></tr>'

    related_sections = []
    for raw_path in related_html_paths:
        path = Path(raw_path)
        title = path.stem.replace("_", " ").replace("-", " ").title()
        rel = _relative_link(path, output_path.parent)
        related_sections.append(
            "<section class=\"embed-card\">"
            f"<h3>{escape(title)}</h3>"
            f"<p><a href=\"{escape(rel)}\">Open {escape(path.name)}</a></p>"
            f"<iframe src=\"{escape(rel)}\" title=\"{escape(title)}\"></iframe>"
            "</section>"
        )
    related_html = "\n".join(related_sections) or (
        '<p class="muted">No additional HTML reports were discovered for this run.</p>'
    )

    network_rel = _relative_link(network_html_path, output_path.parent)
    edges_rel = _relative_link(edges_tsv_path, output_path.parent)
    ecological_edges = evidence_counts.get("co_abundance_guild", 0)
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ContigWeaver Comprehensive Report</title>
  <style>
    :root {{
      --bg: #f4efe6;
      --panel: #fffdf8;
      --ink: #1d2a33;
      --muted: #5b6b74;
      --line: #d8cec0;
      --accent: #146c78;
      --accent-2: #8a5a44;
    }}
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; font-family: Georgia, "Times New Roman", serif; background: linear-gradient(180deg, #f7f1e7 0%, #efe8db 100%); color: var(--ink); }}
    header {{ padding: 48px 24px 28px; background: radial-gradient(circle at top left, rgba(20,108,120,0.18), transparent 34%), linear-gradient(135deg, #18343b, #254b4d); color: #f7f5ef; }}
    header h1 {{ margin: 0 0 10px; font-size: clamp(2rem, 5vw, 3.2rem); line-height: 1.05; }}
    header p {{ margin: 0; max-width: 70ch; color: rgba(247,245,239,0.88); }}
    main {{ padding: 24px; display: grid; gap: 22px; }}
    section {{ background: var(--panel); border: 1px solid var(--line); border-radius: 18px; padding: 20px; box-shadow: 0 18px 40px rgba(29,42,51,0.08); }}
    h2 {{ margin-top: 0; font-size: 1.45rem; }}
    h3 {{ margin-bottom: 8px; }}
    .muted {{ color: var(--muted); }}
    .stats {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 14px; }}
    .stat {{ padding: 16px; border-radius: 14px; background: linear-gradient(180deg, rgba(20,108,120,0.08), rgba(138,90,68,0.08)); border: 1px solid rgba(20,108,120,0.12); }}
    .stat strong {{ display: block; font-size: 1.8rem; margin-bottom: 6px; }}
    .artifact-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 14px; }}
    .artifact {{ border: 1px solid var(--line); border-radius: 14px; padding: 14px; background: #fffaf2; }}
    table {{ width: 100%; border-collapse: collapse; }}
    th, td {{ padding: 10px 12px; text-align: left; border-bottom: 1px solid var(--line); vertical-align: top; }}
    th {{ font-size: 0.9rem; text-transform: uppercase; letter-spacing: 0.04em; color: var(--muted); }}
    a {{ color: var(--accent); }}
    iframe {{ width: 100%; min-height: 640px; border: 1px solid var(--line); border-radius: 14px; background: white; }}
    .embed-grid {{ display: grid; gap: 18px; }}
    .embed-card {{ border-top: 1px solid var(--line); padding-top: 10px; }}
    .callout {{ padding: 14px 16px; border-left: 4px solid var(--accent-2); background: rgba(138,90,68,0.08); border-radius: 10px; }}
    @media (max-width: 700px) {{
      header {{ padding: 32px 18px 20px; }}
      main {{ padding: 16px; }}
      section {{ padding: 16px; }}
      iframe {{ min-height: 420px; }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>ContigWeaver Comprehensive Report</h1>
    <p>Generated {escape(generated_at)}. This report combines the exported network, graph-level evidence summaries, run inputs, and related HTML artifacts so the results are not limited to a graph visualization alone.</p>
  </header>
  <main>
    <section>
      <h2>Run Overview</h2>
      <div class="stats">
        <div class="stat"><strong>{node_count:,}</strong><span>Total nodes</span></div>
        <div class="stat"><strong>{edge_count:,}</strong><span>Total edges</span></div>
        <div class="stat"><strong>{viral_count:,}</strong><span>Viral contigs</span></div>
        <div class="stat"><strong>{ecological_edges:,}</strong><span>Ecological edges</span></div>
      </div>
      <div class="callout">
        Primary outputs: <a href="{escape(network_rel)}">{escape(network_html_path.name)}</a> and <a href="{escape(edges_rel)}">{escape(edges_tsv_path.name)}</a>. This `index.html` is the top-level landing page for the full result set.
      </div>
    </section>

    <section>
      <h2>Inputs Used</h2>
      <table>
        <thead><tr><th>Input</th><th>Path</th><th>Status</th></tr></thead>
        <tbody>
          {input_table}
        </tbody>
      </table>
    </section>

    <section>
      <h2>Evidence Breakdown</h2>
      <div class="artifact-grid">
        <div class="artifact">
          <h3>Edge Types</h3>
          <table>
            <thead><tr><th>Evidence</th><th>Count</th></tr></thead>
            <tbody>{evidence_rows}</tbody>
          </table>
        </div>
        <div class="artifact">
          <h3>Top Connected Contigs</h3>
          <table>
            <thead><tr><th>Contig</th><th>Degree</th></tr></thead>
            <tbody>{degree_rows}</tbody>
          </table>
        </div>
      </div>
    </section>

    <section>
      <h2>Edge Preview</h2>
      <p class="muted">First 15 exported graph edges shown as a quick inspection table.</p>
      <table>
        <thead><tr><th>Source</th><th>Target</th><th>Type</th><th>Weight</th></tr></thead>
        <tbody>{sample_rows}</tbody>
      </table>
    </section>

    <section>
      <h2>Interactive Network</h2>
      <p><a href="{escape(network_rel)}">Open standalone network view</a></p>
      <iframe src="{escape(network_rel)}" title="ContigWeaver interactive network"></iframe>
    </section>

    <section>
      <h2>Related HTML Reports</h2>
      <div class="embed-grid">
        {related_html}
      </div>
    </section>
  </main>
</body>
</html>
"""
    output_path.write_text(html)
    logger.info("Comprehensive HTML report written to %s.", output_path)
    return output_path


def _relative_link(path: Path, base_dir: Path) -> str:
    return relpath(path, start=base_dir).replace("\\", "/")
