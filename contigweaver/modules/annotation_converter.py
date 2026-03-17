import csv
import logging
from collections import defaultdict
from pathlib import Path
import re

logger = logging.getLogger(__name__)
NODE_NUM_RE = re.compile(r"^(NODE_\d+)_")


class ProkkaAnnotationConverter:
    def convert(
        self,
        annotations_input: str | Path,
        output_tsv: str | Path,
        reference_contigs_fasta: str | Path | None = None,
    ) -> Path:
        annotations_input = Path(annotations_input)
        output_tsv = Path(output_tsv)
        reference_map = self._build_reference_id_map(reference_contigs_fasta)

        gff_files = self._discover_gff_files(annotations_input)
        if not gff_files:
            raise ValueError(f"No Prokka GFF files found in: {annotations_input}")

        contig_terms: dict[str, set[str]] = defaultdict(set)
        for gff_path in gff_files:
            tsv_terms = self._parse_prokka_tsv(gff_path.with_suffix(".tsv"))
            self._merge_gff_annotations(gff_path, tsv_terms, contig_terms, reference_map)

        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        with output_tsv.open("w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["Contig_ID", "functional_terms"])
            for contig_id in sorted(contig_terms):
                terms = sorted(contig_terms[contig_id])
                if terms:
                    writer.writerow([contig_id, ",".join(terms)])

        logger.info(
            "Converted Prokka annotations from %d GFF files into %s (%d contigs).",
            len(gff_files),
            output_tsv,
            sum(1 for terms in contig_terms.values() if terms),
        )
        return output_tsv

    def _discover_gff_files(self, annotations_input: Path) -> list[Path]:
        if annotations_input.is_file() and annotations_input.suffix == ".gff":
            return [annotations_input]
        if annotations_input.is_dir():
            return sorted(annotations_input.glob("**/*.gff"))
        return []

    def _parse_prokka_tsv(self, tsv_path: Path) -> dict[str, set[str]]:
        if not tsv_path.exists():
            return {}

        terms_by_locus: dict[str, set[str]] = defaultdict(set)
        with tsv_path.open(newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                locus_tag = (row.get("locus_tag") or "").strip()
                if not locus_tag:
                    continue

                for field in ("gene", "EC_number", "COG", "product"):
                    value = (row.get(field) or "").strip()
                    if value:
                        terms_by_locus[locus_tag].update(self._normalize_terms(value))

        return terms_by_locus

    def _merge_gff_annotations(
        self,
        gff_path: Path,
        tsv_terms: dict[str, set[str]],
        contig_terms: dict[str, set[str]],
        reference_map: dict[str, str],
    ) -> None:
        with gff_path.open() as fh:
            for raw_line in fh:
                if raw_line.startswith("#"):
                    continue

                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "CDS":
                    continue

                contig_id = self._canonicalize_contig_id(fields[0], reference_map)
                attrs = self._parse_attributes(fields[8])
                locus_tag = attrs.get("locus_tag") or attrs.get("ID")

                terms: set[str] = set()
                for key in ("gene", "Name", "product", "eC_number"):
                    value = attrs.get(key)
                    if value:
                        terms.update(self._normalize_terms(value))

                for value in attrs.get("db_xref", "").split(","):
                    value = value.strip()
                    if not value:
                        continue
                    if ":" in value:
                        _, value = value.split(":", 1)
                    terms.update(self._normalize_terms(value))

                if locus_tag and locus_tag in tsv_terms:
                    terms.update(tsv_terms[locus_tag])

                contig_terms[contig_id].update(terms)

    @staticmethod
    def _parse_attributes(raw_attrs: str) -> dict[str, str]:
        attrs: dict[str, str] = {}
        for item in raw_attrs.split(";"):
            if "=" not in item:
                continue
            key, value = item.split("=", 1)
            attrs[key] = value
        return attrs

    def _normalize_terms(self, raw_value: str) -> set[str]:
        raw_value = raw_value.strip()
        if not raw_value or raw_value.lower() == "hypothetical protein":
            return set()

        terms = {raw_value.upper()}
        expanded = {raw_value}
        for sep in ("/", "|", ";"):
            next_values: set[str] = set()
            for value in expanded:
                next_values.update(part.strip() for part in value.split(sep) if part.strip())
            expanded = next_values or expanded

        for value in expanded:
            upper = value.upper()
            terms.add(upper)
            terms.update(self._keyword_terms(upper))

        return {term for term in terms if term}

    @staticmethod
    def _keyword_terms(value: str) -> set[str]:
        keyword_map = {
            "GLYCOLYSIS": "GLYCOLYSIS",
            "CITRATE SYNTHASE": "TCA",
            "ISOCITRATE DEHYDROGENASE": "TCA",
            "MALATE DEHYDROGENASE": "TCA",
            "SUCCINATE DEHYDROGENASE": "TCA",
            "NITRATE REDUCTASE": "NITRATE-DEG",
            "NITRITE REDUCTASE": "NITRATE-DEG",
            "NITRIC OXIDE REDUCTASE": "NITRATE-DEG",
            "NITROUS-OXIDE REDUCTASE": "NITRATE-DEG",
            "SULFATE REDUCTASE": "SULFATE-REDUCTION",
            "SULFITE REDUCTASE": "SULFATE-REDUCTION",
            "ADENYLYLSULFATE REDUCTASE": "SULFATE-REDUCTION",
            "SULFATE ADENYLYLTRANSFERASE": "SULFATE-REDUCTION",
            "PHOTOSYSTEM": "PHOTOSYSTEM",
            "RIBULOSE-BISPHOSPHATE CARBOXYLASE": "PHOTOSYSTEM",
        }
        return {term for needle, term in keyword_map.items() if needle in value}

    @staticmethod
    def _build_reference_id_map(
        reference_contigs_fasta: str | Path | None,
    ) -> dict[str, str]:
        if reference_contigs_fasta is None:
            return {}

        reference_map: dict[str, str] = {}
        with Path(reference_contigs_fasta).open() as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                contig_id = line[1:].strip().split()[0]
                match = NODE_NUM_RE.match(contig_id)
                if match:
                    reference_map[match.group(1)] = contig_id
        return reference_map

    @staticmethod
    def _canonicalize_contig_id(contig_id: str, reference_map: dict[str, str]) -> str:
        match = NODE_NUM_RE.match(contig_id)
        if not match:
            return contig_id
        return reference_map.get(match.group(1), contig_id)


def prepare_annotations_input(
    annotations_input: str | Path | None,
    work_dir: str | Path,
    reference_contigs_fasta: str | Path | None = None,
) -> Path | None:
    if annotations_input is None:
        return None

    annotations_input = Path(annotations_input)
    if annotations_input.is_file():
        return annotations_input
    if annotations_input.is_dir():
        output_tsv = Path(work_dir) / "converted_annotations.tsv"
        return ProkkaAnnotationConverter().convert(
            annotations_input,
            output_tsv,
            reference_contigs_fasta=reference_contigs_fasta,
        )
    raise FileNotFoundError(f"Annotations input not found: {annotations_input}")
