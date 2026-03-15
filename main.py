#!/usr/bin/env python3
"""
ContigWeaver — top-level entry point.

Run the pipeline directly without installing the package:

    python main.py --gfa assembly.gfa.gz --contigs contigs.fa.gz \\
                   --viral-contigs viral.fa --output-dir out/

This is a thin shim; all logic lives in contigweaver/pipeline.py.
"""
import sys
from pathlib import Path

# Allow running from the project root without `pip install -e .`
sys.path.insert(0, str(Path(__file__).parent))

from contigweaver.pipeline import main

if __name__ == "__main__":
    sys.exit(main())
