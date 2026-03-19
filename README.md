# ContigWeaver

> **Weave inter-contig interaction networks from low-depth metagenomic data — no MAG binning required.**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## Overview

Most metagenomic binning workflows collapse complex microbial communities into MAGs (Metagenome-Assembled Genomes). This works well for high-coverage datasets, but **low-depth marine and environmental metagenomes** often yield fragmented assemblies where bins are incomplete, chimeric, or simply absent.

**ContigWeaver** takes a different approach: instead of binning, it builds a **multi-evidence interaction network** directly from the assembly graph, CRISPR spacer targeting, co-abundance patterns, and metabolic complementarity. Physical evidence starts from native GFA segment IDs, while contig-level evidence (`NODE_*`) is bridged back into the graph with sequence-derived `segment_membership` anchors so the layers stay connected on real SPAdes assemblies.

The result is a browsable interactive HTML network — a "contig web" that reveals phage–host linkages, co-occurring community members, and metabolic guilds even when coverage is too low to bin.

---

## Key Features

| Feature | Description |
|---|---|
| **Binning-free** | Works directly on contigs from any assembler (SPAdes, MEGAHIT) |
| **Bridge-aware integration** | Reconciles `NODE_*` contigs back to numeric GFA segments with exact or substring sequence anchors |
| **Seven evidence layers** | Physical overlap · Segment membership · CRISPR targeting · Co-abundance · Taxonomic match · Functional operon · Bin membership |
| **Gzip-transparent** | Reads `.gfa.gz` and `.fasta.gz` directly — no manual decompression |
| **Scalable HTML output** | Full TSV export plus an automatically focused HTML view for very large graphs |
| **Modular pipeline** | Stage 1 runs with just GFA + contigs; Stage 2 accepts either a ready TSV or a Prokka directory |
| **Low dependencies** | Pure Python: networkx, pyvis, pandas, scipy |

---

## Evidence Layers

### Stage 1 — Physical Evidence (always required)

| Module | Source | Edge type |
|---|---|---|
| **GFA Parser** | SPAdes / MEGAHIT assembly graph (`.gfa`) | `physical_overlap` |
| **Contig Reconciler** | Sequence anchors between FASTA contigs and GFA segments | `segment_membership` |
| **CRISPR-Phage Miner** | MinCED spacer prediction + BLAST vs. viral contigs | `crispr_targeting` |

### Stage 2 — Ecological + Annotation + Binning Evidence (optional)

| Module | Source | Edge type |
|---|---|---|
| **Co-Abundance Correlator** | Per-sample coverage table (Spearman ρ ≥ 0.85) | `co_abundance_guild` |
| **Pathway Complementarity** | Functional annotations (TSV or converted Prokka directory) | attribute on co-abundance edges |
| **Annotation Miner** | Taxonomy + functional TSV (`--annotation-data`) | `taxonomic_match`, `functional_operon` |
| **Binning Miner** | Contig-to-bin TSV (`--binning`) | `bin_membership` + rescued bin metadata on nodes |

---

## Installation

```bash
git clone https://github.com/your-org/contigweaver.git
cd contigweaver
pip install -e .
```

Or run it directly without installing:

```bash
python main.py --help
```

**Runtime dependencies** (automatically installed):

```
networkx>=3.0
pyvis>=0.3.2
pandas>=2.0
scipy>=1.11
```

**External tools** (optional — only for CRISPR evidence):

```bash
# MinCED — CRISPR array finder
conda install -c bioconda minced

# BLAST+ — spacer-vs-virus alignment
conda install -c bioconda blast
```

---

## Quick Start

### Stage 1 only (GFA + CRISPR)

```bash
python main.py \
  --gfa "input/SPAdes_Asm/SPAdes-sponge_brin.assembly.gfa.gz" \
  --contigs "input/SPAdes_Asm/SPAdes-sponge_brin.contigs.fa" \
  --viral-contigs "input/SPAdes_Genomad/Galaxy48-[geNomad on dataset 43_ virus fasta].fasta" \
  --output-dir runs/
```

### Stage 1 + Stage 2 (full pipeline)

```bash
python main.py \
  --gfa "input/SPAdes_Asm/SPAdes-sponge_brin.assembly.gfa.gz" \
  --contigs "input/SPAdes_Asm/SPAdes-sponge_brin.contigs.fa" \
  --viral-contigs "input/SPAdes_Genomad/Galaxy48-[geNomad on dataset 43_ virus fasta].fasta" \
  --coverage "tmp_real_annotation_eval/coverage_from_headers.tsv" \
  --annotations "input/SPAdes_Prokka" \
  --annotation-data "annotation_data.tsv" \
  --binning "binning.tsv" \
  --output-dir runs/ \
  --verbose
```

### Comprehensive command example

```bash
python main.py \
  --gfa "input/SPAdes_Asm/SPAdes-sponge_brin.assembly.gfa.gz" \
  --contigs "input/SPAdes_Asm/SPAdes-sponge_brin.contigs.fa" \
  --viral-contigs "input/SPAdes_Genomad/Galaxy48-[geNomad on dataset 43_ virus fasta].fasta" \
  --coverage "tmp_real_annotation_eval/coverage_from_headers.tsv" \
  --annotations "input/SPAdes_Prokka" \
  --annotation-data "annotation_data.tsv" \
  --binning "binning.tsv" \
  --minced-bin minced \
  --blastn-bin blastn \
  --makeblastdb-bin makeblastdb \
  --spearman-threshold 0.85 \
  --p-value-threshold 0.01 \
  --output-dir runs/ \
  --verbose
```

Each execution now creates an isolated run folder under `runs/`, for example:

- `runs/run_20260319_082502_stage1_a1b2c3d4/`
- `runs/run_20260319_083044_stage2_e5f6g7h8/`

Use `runs/latest/` to jump to the newest run, then open `runs/latest/contigweaver_network.html` in any browser to explore the network.

For large real assemblies, ContigWeaver always writes the **full TSV** and then writes a **focused HTML subgraph** around non-physical evidence so browser rendering stays practical.

---

## Input Files

| Flag | Format | Description |
|---|---|---|
| `--gfa` | GFA v1 (plain or `.gz`) | Assembly graph from SPAdes or MEGAHIT |
| `--contigs` | FASTA (plain or `.gz`) | All assembled contigs |
| `--viral-contigs` | FASTA | Viral / phage contig subset (see [Identifying viral contigs](#identifying-viral-contigs)) |
| `--coverage` *(optional)* | TSV: `Contig_ID \| sample_1 \| sample_2 \| …` | Per-sample coverage — minimum 3 samples required |
| `--annotations` *(optional)* | TSV or Prokka directory | TSV: `Contig_ID \| KO_terms \| MetaCyc_terms` or `Contig_ID \| functional_terms`; directories of Prokka `.gff`/`.tsv` files are converted automatically |
| `--annotation-data` *(optional)* | TSV | Annotation Miner input: `Contig_ID` plus optional taxonomy/functional/CRISPR columns |
| `--binning` *(optional)* | TSV | Binning Miner input: `Contig_ID \| Bin_ID` |

### Annotation Miner TSV format (`--annotation-data`)

- Required: `Contig_ID`
- Optional taxonomy columns: `taxonomy_label`, `taxonomy_rank`, `taxonomy_confidence`
- Optional functional columns: `functional_terms`, `cas_genes`
- Optional CRISPR columns: `has_crispr`, `spacer_count`

Example:

```tsv
Contig_ID	taxonomy_label	taxonomy_rank	taxonomy_confidence	functional_terms	cas_genes	has_crispr	spacer_count
NODE_1	Bacillus	genus	0.98	CRISPR repeat,defense		true	3
NODE_2	Bacillus	genus	0.97	DNA repair	Cas9	false	0
```

### Binning TSV format (`--binning`)

- Required: `Contig_ID`, `Bin_ID`

Example:

```tsv
Contig_ID	Bin_ID
NODE_1	bin_001
NODE_2	bin_001
NODE_9	bin_004
```

---

## Output Files

| File | Description |
|---|---|
| `contigweaver_edges.tsv` | Full edge list with Source, Target, Evidence_Type, weight, attributes |
| `contigweaver_network.html` | **Interactive network** — full graph if small; otherwise a focused ego-view around the top actor with evidence-layer controls and auto summary |
| `index.html` | Comprehensive run report with overview, evidence breakdown, edge preview, network embed, related HTML reports, and automated multi-evidence Key Findings |
| `run.log` | Per-run log file with start/end time, arguments, module progress, stats, warnings, and errors |
| `spacers.fasta` | Extracted CRISPR spacers (Stage 1) |
| `spacer_vs_viral.tsv` | BLAST hits: spacer → viral contig (identity ≥ 95 %, coverage ≥ 90 %) |
| `minced_output.txt` | Raw MinCED predictions |
| `workdir/converted_annotations.tsv` | Generated automatically when `--annotations` points to a Prokka directory |

All files above are written inside the run-specific folder (`runs/run_*`) and never in the repository root.

### Large-graph behavior

- `contigweaver_edges.tsv` always contains the complete graph.
- `contigweaver_network.html` switches to a focused view when the graph exceeds the HTML budget.
- Physical hairball suppression is applied in HTML mode: `physical_overlap` edges are retained only when they connect viral and non-viral nodes.
- Current default HTML budget: `1500` nodes and `4000` edges.
- Focus center is selected from the highest-degree viral node (fallback: global top-degree node), then expanded using an ego-network (`radius=2`).
- Default edge-layer toggles in HTML: `crispr_targeting`, `segment_membership`, `co_abundance_guild` ON; `physical_overlap` OFF.
- Focus priority: non-physical evidence first (`segment_membership`, `crispr_targeting`, `co_abundance_guild`, `functional_operon`, `taxonomic_match`), then local graph neighborhood.

### Annotation directory support

- `--annotations` can point directly to a Prokka results directory such as `input/SPAdes_Prokka/`.
- Contig IDs from per-bin Prokka files are canonicalized back to the Stage 1 assembly FASTA by shared `NODE_<index>`.
- Prokka `gene`, `product`, `EC_number`, and `COG` terms are aggregated per contig into `functional_terms`.
- The generated TSV is written under `workdir/converted_annotations.tsv` and then consumed by Stage 2 exactly like a hand-written annotations table.

---

## Identifying Viral Contigs

ContigWeaver needs a FASTA of putative viral contigs for CRISPR targeting evidence. Options:

### Option A — geNomad (recommended)

```bash
conda install -c conda-forge -c bioconda genomad
genomad download-database ./genomad_db
genomad end-to-end \
  results/Assembly/SPAdes/SPAdes-sponge_brin.contigs.fa.gz \
  genomad_output/ genomad_db/
# result: genomad_output/*/sponge_brin_contigs_virus.fna
```

### Option B — VirSorter2

```bash
conda install -c conda-forge -c bioconda virsorter2
virsorter run -w vs2_output \
  -i results/Assembly/SPAdes/SPAdes-sponge_brin.contigs.fa.gz \
  --min-length 1500 -j 8 all
# result: vs2_output/final-viral-combined.fa
```

### Option C — Heuristic (no install, quick)

```bash
python3 - <<'EOF'
import gzip, re
in_fa  = "results/Assembly/SPAdes/SPAdes-sponge_brin.contigs.fa.gz"
out_fa = "viral_contigs_heuristic.fa"
MIN_LEN, MAX_LEN, MIN_COV = 1500, 50_000, 100.0
written = 0
with gzip.open(in_fa, "rt") as fh, open(out_fa, "w") as out:
    header, seq = "", []
    for line in fh:
        if line.startswith(">"):
            if header:
                m = re.search(r"length_(\d+)_cov_([\d.]+)", header)
                if m and MIN_LEN <= int(m.group(1)) <= MAX_LEN and float(m.group(2)) >= MIN_COV:
                    out.write(header + "".join(seq) + "\n"); written += 1
            header, seq = line, []
        else:
            seq.append(line)
print(f"{written} putative viral contigs → {out_fa}")
EOF
```

---

## CLI Reference

```
usage: contigweaver [-h] --gfa FILE --contigs FILE --viral-contigs FILE
                    [--coverage FILE] [--annotations FILE]
                    [--annotation-data FILE] [--binning FILE]
                    [--output-dir DIR]
                    [--minced-bin PATH] [--blastn-bin PATH] [--makeblastdb-bin PATH]
                    [--spearman-threshold FLOAT] [--p-value-threshold FLOAT]
                    [-v]

options:
  --gfa FILE                 Assembly graph (GFA v1, plain or .gz)
  --contigs FILE             All contigs FASTA (plain or .gz)
  --viral-contigs FILE       Viral / phage contigs FASTA
  --coverage FILE            Per-sample coverage table TSV  [Stage 2]
  --annotations FILE         Functional annotations TSV or Prokka directory  [Stage 2]
  --annotation-data FILE     Annotation Miner TSV  [Stage 2]
  --binning FILE             Binning assignments TSV with Contig_ID/Bin_ID  [Stage 2]
  --output-dir DIR           Root output directory for run_* folders (default: runs)
  --minced-bin PATH          Path to minced binary (default: minced)
  --blastn-bin PATH          Path to blastn binary (default: blastn)
  --makeblastdb-bin PATH     Path to makeblastdb binary (default: makeblastdb)
  --spearman-threshold FLOAT Minimum |ρ| for co-abundance edge (default: 0.85)
  --p-value-threshold FLOAT  Maximum p-value for co-abundance edge (default: 0.01)
  -v, --verbose              Enable debug logging
```

---

## Python API

```python
from contigweaver.pipeline import ContigWeaverPipeline

pipeline = ContigWeaverPipeline(output_dir="my_output/")

# Stage 1
pipeline.run_stage1(
    gfa_path="assembly.gfa.gz",
    contigs_fasta="contigs.fa.gz",
    viral_contigs_fasta="viral.fa",
)

# Stage 2 (optional)
pipeline.run_stage2(
    coverage_tsv="coverage.tsv",
    annotations_tsv="annotations.tsv",   # optional TSV or Prokka directory
    annotation_data_tsv="annotation_data.tsv",  # optional taxonomy/functional layer
    binning_tsv="binning.tsv",  # optional contig-to-bin mapping + rescue
)

# Access the networkx graph directly
graph = pipeline.graph
print(f"{graph.number_of_nodes()} contigs, {graph.number_of_edges()} interactions")
```

---

## Project Structure

```
contigweaver/
├── __init__.py
├── __main__.py
├── pipeline.py                  # ContigWeaverPipeline + CLI entry point
├── modules/
│   ├── gfa_parser.py            # Module 1 — GFA physical overlap
│   ├── crispr_miner.py          # Module 2 — CRISPR-phage targeting
│   ├── graph_exporter.py        # Module 3 — TSV + HTML export
│   ├── ecological_miner.py      # Modules 4–6 — co-abundance + metabolic
│   ├── annotation_miner.py      # Stage 2A — taxonomy + functional operon links
│   └── binning_miner.py         # Stage 2B — bin integration + graph-aware rescue
├── tests/
│   ├── conftest.py              # Shared fixtures (mini SPAdes-format GFA)
│   ├── test_gfa_parser.py
│   ├── test_crispr_miner.py
│   ├── test_graph_exporter.py
│   ├── test_ecological_miner.py
│   ├── test_annotation_miner.py
│   ├── test_binning_miner.py
│   └── test_pipeline_integration.py
└── data/
    └── example/                 # Small synthetic test dataset
```

---

## Running Tests

```bash
pip install -e ".[dev]"
pytest contigweaver/tests/ -v
```

The integration test suite includes a `real_data` smoke test that reads the actual SPAdes GFA from `results/` (auto-skipped if not present):

```
contigweaver/tests/test_pipeline_integration.py::TestRealDataSmoke::test_real_spades_gfa_parses
```

---

## Background

This tool was developed for the analysis of **low-depth marine sponge metagenomes** (sample: *sponge_brin*) where conventional binning tools (CONCOCT, MaxBin2, MetaBAT2) produce fragmented, incomplete MAGs. By reasoning at the contig level and weaving together multiple lines of evidence, ContigWeaver recovers ecological signals — particularly phage–host and metabolic guild relationships — that are invisible to binning-first pipelines.

---

## License

MIT © ContigWeaver Contributors
