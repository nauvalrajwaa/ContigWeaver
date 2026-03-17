from pathlib import Path

import networkx as nx

from contigweaver.modules.annotation_converter import (
    ProkkaAnnotationConverter,
    prepare_annotations_input,
)
from contigweaver.modules.ecological_miner import PathwayComplementarityChecker
from contigweaver.pipeline import ContigWeaverPipeline


def _write_mock_prokka_dir(base_dir: Path) -> Path:
    prokka_dir = base_dir / "mock_prokka" / "bin_1"
    prokka_dir.mkdir(parents=True)

    (prokka_dir / "bin_1.gff").write_text(
        "##gff-version 3\n"
        "NODE_1\tProkka\tCDS\t1\t900\t.\t+\t0\tID=LOC_0001;locus_tag=LOC_0001;gene=gltA;product=citrate synthase;db_xref=COG:COG0372\n"
        "NODE_1\tProkka\tCDS\t950\t1500\t.\t+\t0\tID=LOC_0002;locus_tag=LOC_0002;product=hypothetical protein\n"
        "NODE_2\tProkka\tCDS\t1\t1200\t.\t-\t0\tID=LOC_0003;locus_tag=LOC_0003;gene=mdh;product=malate dehydrogenase;eC_number=1.1.1.37\n"
    )
    (prokka_dir / "bin_1.tsv").write_text(
        "locus_tag\tftype\tlength_bp\tgene\tEC_number\tCOG\tproduct\n"
        "LOC_0001\tCDS\t900\tgltA\t2.3.3.1\tCOG0372\tcitrate synthase\n"
        "LOC_0003\tCDS\t1200\tmdh\t1.1.1.37\tCOG0173\tmalate dehydrogenase\n"
    )

    return prokka_dir.parent


def test_prokka_converter_writes_contig_level_functional_terms(tmp_path: Path):
    prokka_root = _write_mock_prokka_dir(tmp_path)
    output_tsv = tmp_path / "converted.tsv"

    result = ProkkaAnnotationConverter().convert(prokka_root, output_tsv)

    assert result == output_tsv
    content = output_tsv.read_text()
    assert "Contig_ID\tfunctional_terms" in content
    assert "NODE_1" in content
    assert "NODE_2" in content
    assert "TCA" in content
    assert "TCA" in content


def test_prepare_annotations_input_accepts_directory(tmp_path: Path):
    prokka_root = _write_mock_prokka_dir(tmp_path)

    converted = prepare_annotations_input(prokka_root, tmp_path / "workdir")

    assert converted is not None
    assert converted.exists()
    checker = PathwayComplementarityChecker()
    checker.load_annotations(converted)
    assert checker.check_pair("NODE_1", "NODE_1") is True


def test_pipeline_stage2_accepts_prokka_directory(tmp_path: Path):
    prokka_root = _write_mock_prokka_dir(tmp_path)
    coverage_tsv = tmp_path / "coverage.tsv"
    coverage_tsv.write_text(
        "Contig_ID\ts1\ts2\ts3\n"
        "NODE_1\t10\t20\t30\n"
        "NODE_2\t11\t21\t31\n"
        "NODE_3\t1\t5\t2\n"
    )

    pipeline = ContigWeaverPipeline(output_dir=tmp_path / "out")
    pipeline.graph.add_node("NODE_1", node_type="unknown", length=1000)
    pipeline.graph.add_node("NODE_2", node_type="unknown", length=1200)

    pipeline.run_stage2(coverage_tsv=coverage_tsv, annotations_tsv=prokka_root)

    converted = tmp_path / "out" / "workdir" / "converted_annotations.tsv"
    assert converted.exists()
    eco_edges = [
        data
        for _, _, data in pipeline.graph.edges(data=True)
        if data.get("type") == "co_abundance_guild"
    ]
    assert eco_edges
    assert eco_edges[0]["metabolic_match"] is True
