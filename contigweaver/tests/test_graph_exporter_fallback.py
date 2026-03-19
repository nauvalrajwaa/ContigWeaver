import builtins

import networkx as nx

from contigweaver.modules.graph_exporter import GraphExporter


def test_export_html_writes_fallback_when_pyvis_missing(tmp_path, monkeypatch):
    graph = nx.MultiGraph()
    graph.add_node("host_1", node_type="unknown", length=1200)
    graph.add_node("virus_1", node_type="viral", length=700)
    graph.add_edge("host_1", "virus_1", type="crispr_targeting", identity=98.2, coverage=0.97)

    original_import = builtins.__import__

    def patched_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "pyvis.network" or name == "pyvis":
            raise ImportError("simulated missing pyvis")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", patched_import)

    output_html = tmp_path / "network.html"
    GraphExporter(graph).export_html(output_html, viral_contigs={"virus_1"})

    content = output_html.read_text()
    assert output_html.exists()
    assert "Fallback HTML" in content
    assert "pyvis" in content
    assert "Evidence Breakdown" in content
