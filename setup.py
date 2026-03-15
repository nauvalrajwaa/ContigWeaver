"""
setup.py for ContigWeaver — replaces pyproject.toml for plain-Python packaging.

Install (editable):
    pip install -e .

Install (with dev tools):
    pip install -e ".[dev]"
"""
from setuptools import setup, find_packages

setup(
    name="contigweaver",
    version="1.0.0",
    description=(
        "Reconstruct inter-contig interaction networks from low-depth "
        "metagenomic data without MAG binning"
    ),
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.10",
    keywords=[
        "metagenomics",
        "bioinformatics",
        "assembly graph",
        "CRISPR",
        "network",
        "contig weaving",
        "binning-free",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: MIT License",
    ],
    packages=find_packages(include=["contigweaver", "contigweaver.*"]),
    install_requires=[
        "networkx>=3.0",
        "pyvis>=0.3.2",
        "pandas>=2.0",
        "scipy>=1.11",
    ],
    extras_require={
        "dev": [
            "pytest>=7.4",
            "pytest-cov>=4.1",
        ],
    },
    entry_points={
        "console_scripts": [
            "contigweaver=contigweaver.pipeline:main",
        ],
    },
)
