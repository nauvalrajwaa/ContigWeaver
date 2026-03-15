"""
Shared pytest fixtures for ContigNexus integration testing.

Provides a compact, self-consistent fixture set that mirrors the real SPAdes
output format (integer segment IDs, DP:f: depth tags, 55M overlaps, .gz files)
found in results/Assembly/SPAdes/SPAdes-sponge_brin.assembly.gfa.gz.
"""
import gzip
import io
import textwrap
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Segment sequences sampled from the first 20 S-lines of the real SPAdes GFA
# (results/Assembly/SPAdes/SPAdes-sponge_brin.assembly.gfa.gz).
# IDs 5, 7, 11 are bacterial contigs; IDs 19, 37 are high-depth (viral proxy).
# ---------------------------------------------------------------------------

SPADES_MINI_GFA = textwrap.dedent("""\
    H\tVN:Z:1.2\tsp:Z:SPAdes-4.1.0
    S\t5\tTTATTTGATGCTAAAAGCGCCGAATTATTCACGTTTTTCTGATCGAATTTTTTCCAGCCCTCTAGCGTTCTTTTTCGGAACTAATCGCGATCGCGCAATGTGCCCCTGGTGTTAATTTTTGGAACTACGAAAAAGGAGAAAATAAGCGGCAAAATGCGAGAATTCCCTGCAAACC\tDP:f:1.91195\tKC:i:304
    S\t7\tACGAGAGGTGGGATTGGTGGGATGTCAAGTGGGATTGGTGGGACGAGAGGTGGGACTGGTGGAATGAGGAGTGGGATTGGTGGGACGAGAGGTGGGATTGGTGGGATGA\tDP:f:63.7593\tKC:i:3443
    S\t11\tGAACGAAGGGAAATAATACCGCAGTATTTTTTTCACACtCAACTACCTTTAGGTTCTGATCATTCCCATCAACATCAAATTCACAAAT\tDP:f:1.18038\tKC:i:373
    S\t17\tTTCACCACTAATTGAGAATAAATAAATTGTAGAATTATCATTGTGAAACATTCAGTATATTAATTTTTTTTAATCCAAAATGATAATTCGAAAAATAAGGTACCTGTC\tDP:f:1.36607\tKC:i:765
    S\t19\tGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGGGTGGGTGTGTGTGGGTGTGTGTGTGTGTGTGTGTGTGTGTG\tDP:f:12.7143\tKC:i:356
    S\t37\tGGCCTTTCCTGACATCATGTCTCAACAGAGAATTTCAGAAGAGGCAAGTAGTACCTCCTCACCCCCTTCTCCTCGTCCTCCTTCTCCTCCTCGTCCTCCTCGTCCTCCTCCTCCTCACATCC\tDP:f:14.9354\tKC:i:16638
    L\t5\t+\t11\t+\t55M
    L\t11\t+\t17\t-\t55M
    L\t5\t-\t19\t+\t55M
    L\t7\t+\t37\t+\t55M
    L\t19\t+\t37\t-\t55M
""")

# Mini contigs FASTA (SPAdes NODE_<n>_length_<bp>_cov_<depth> naming)
MINI_CONTIGS_FASTA = textwrap.dedent("""\
    >NODE_5_length_172_cov_1.91195
    TTATTTGATGCTAAAAGCGCCGAATTATTCACGTTTTTCTGATCGAATTTTTTCCAGCCCTCTAGCGTTCTTTTTCGGAACT
    AATCGCGATCGCGCAATGTGCCCCTGGTGTTAATTTTTGGAACTACGAAAAAGGAGAAAATAAGCGGCAAAATGCGAGAATT
    CCCTGCAAACC
    >NODE_7_length_110_cov_63.7593
    ACGAGAGGTGGGATTGGTGGGATGTCAAGTGGGATTGGTGGGACGAGAGGTGGGACTGGTGGAATGAGGAGTGGGATTGGTG
    GGACGAGAGGTGGGATTGGTGGGATGA
    >NODE_11_length_88_cov_1.18038
    GAACGAAGGGAAATAATACCGCAGTATTTTTTTCACACTCAACTACCTTTAGGTTCTGATCATTCCCATCAACATCAAATTCA
    CAAAAT
    >NODE_17_length_108_cov_1.36607
    TTCACCACTAATTGAGAATAAATAAATTGTAGAATTATCATTGTGAAACATTCAGTATATTAATTTTTTTTAATCCAAAATGA
    TAATTCGAAAAATAAGGTACCTGTC
    >NODE_19_length_82_cov_12.7143
    GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGGGTGGGTGTGTGTGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
    >NODE_37_length_124_cov_14.9354
    GGCCTTTCCTGACATCATGTCTCAACAGAGAATTTCAGAAGAGGCAAGTAGTACCTCCTCACCCCCTTCTCCTCGTCCTCCTT
    CTCCTCCTCGTCCTCCTCGTCCTCCTCCTCCTCACATCC
""")

# Viral contigs: segments 7 (DP=63.76) and 37 (DP=14.93) are high-depth, used as viral proxy
MINI_VIRAL_FASTA = textwrap.dedent("""\
    >NODE_7_length_110_cov_63.7593
    ACGAGAGGTGGGATTGGTGGGATGTCAAGTGGGATTGGTGGGACGAGAGGTGGGACTGGTGGAATGAGGAGTGGGATTGGTG
    GGACGAGAGGTGGGATTGGTGGGATGA
    >NODE_37_length_124_cov_14.9354
    GGCCTTTCCTGACATCATGTCTCAACAGAGAATTTCAGAAGAGGCAAGTAGTACCTCCTCACCCCCTTCTCCTCGTCCTCCTT
    CTCCTCCTCGTCCTCCTCGTCCTCCTCCTCCTCACATCC
""")

# Minimal coverage table (3 samples minimum)
MINI_COVERAGE_TSV = textwrap.dedent("""\
    Contig_ID\tsample_1\tsample_2\tsample_3
    NODE_5_length_172_cov_1.91195\t1.9\t2.1\t1.8
    NODE_7_length_110_cov_63.7593\t63.0\t65.2\t61.8
    NODE_11_length_88_cov_1.18038\t1.2\t1.1\t1.3
    NODE_17_length_108_cov_1.36607\t1.4\t1.3\t1.5
    NODE_19_length_82_cov_12.7143\t12.7\t12.8\t12.6
    NODE_37_length_124_cov_14.9354\t14.9\t15.1\t14.7
""")

# Minimal functional annotations (KO + MetaCyc terms)
MINI_ANNOTATIONS_TSV = textwrap.dedent("""\
    Contig_ID\tKO_terms\tMetaCyc_terms
    NODE_5_length_172_cov_1.91195\tK00001;K00002\tPWY-5345
    NODE_7_length_110_cov_63.7593\tK00003\t
    NODE_11_length_88_cov_1.18038\tK00002;K00004\tPWY-5345;PWY-6737
    NODE_17_length_108_cov_1.36607\t\t
    NODE_19_length_82_cov_12.7143\tK00005\tPWY-7229
    NODE_37_length_124_cov_14.9354\tK00001\tPWY-5345
""")

# ------------------------------------------------------------------
# Fixtures
# ------------------------------------------------------------------

@pytest.fixture
def mini_gfa_file(tmp_path: Path) -> Path:
    """Plain-text GFA file using real SPAdes integer-ID format."""
    p = tmp_path / "mini_spades.gfa"
    p.write_text(MINI_CONTIGS_FASTA.replace(MINI_CONTIGS_FASTA, SPADES_MINI_GFA))
    p.write_text(SPADES_MINI_GFA)
    return p


@pytest.fixture
def mini_gfa_gz_file(tmp_path: Path) -> Path:
    """Gzip-compressed GFA file — mirrors the real results/ file."""
    p = tmp_path / "mini_spades.gfa.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(SPADES_MINI_GFA)
    return p


@pytest.fixture
def mini_contigs_fasta(tmp_path: Path) -> Path:
    p = tmp_path / "contigs.fa"
    p.write_text(MINI_CONTIGS_FASTA)
    return p


@pytest.fixture
def mini_viral_fasta(tmp_path: Path) -> Path:
    p = tmp_path / "viral_contigs.fa"
    p.write_text(MINI_VIRAL_FASTA)
    return p


@pytest.fixture
def mini_coverage_tsv(tmp_path: Path) -> Path:
    p = tmp_path / "coverage.tsv"
    p.write_text(MINI_COVERAGE_TSV)
    return p


@pytest.fixture
def mini_annotations_tsv(tmp_path: Path) -> Path:
    p = tmp_path / "annotations.tsv"
    p.write_text(MINI_ANNOTATIONS_TSV)
    return p


# ------------------------------------------------------------------
# Marker for tests that need the real results/ directory
# ------------------------------------------------------------------

REAL_GFA = Path(__file__).parent.parent.parent / "results/Assembly/SPAdes/SPAdes-sponge_brin.assembly.gfa.gz"

real_data = pytest.mark.skipif(
    not REAL_GFA.exists(),
    reason="Real results/ data not available",
)
