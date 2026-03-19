import csv
import logging
from collections import defaultdict
from pathlib import Path
import re

logger = logging.getLogger(__name__)
NODE_NUM_RE = re.compile(r"^(NODE_\d+)_")

METHOD_ALIASES: dict[str, tuple[str, ...]] = {
    "dastool": ("dastool",),
    "metabat2": ("metabat2", "metabat"),
    "maxbin2": ("maxbin2", "maxbin"),
    "concoct": ("concoct",),
}


def _normalize_for_match(value: str) -> str:
    return "".join(ch for ch in value.lower() if ch.isalnum())


def _detect_method_from_text(value: str) -> str | None:
    normalized = _normalize_for_match(value)
    for method, aliases in METHOD_ALIASES.items():
        if any(alias in normalized for alias in aliases):
            return method
    return None


def detect_method_from_text(value: str) -> str | None:
    return _detect_method_from_text(value)


def detect_methods_from_inputs(
    binning_input: str | Path | list[str | Path] | tuple[str | Path, ...] | None,
) -> set[str]:
    if binning_input is None:
        return set()

    if isinstance(binning_input, (str, Path)):
        values = [str(binning_input)]
    else:
        values = [str(value) for value in binning_input]

    methods = {
        method
        for value in values
        if (method := _detect_method_from_text(value)) is not None
    }
    return methods


class ProkkaAnnotationConverter:
    def convert(
        self,
        annotations_input: str | Path,
        output_tsv: str | Path,
        reference_contigs_fasta: str | Path | None = None,
        binning_input: str | Path | list[str | Path] | tuple[str | Path, ...] | None = None,
        method_selection: set[str] | None = None,
    ) -> Path:
        annotations_input = Path(annotations_input)
        output_tsv = Path(output_tsv)
        reference_map = self._build_reference_id_map(reference_contigs_fasta)

        gff_files = self._discover_gff_files(annotations_input)
        if not gff_files:
            raise ValueError(f"No Prokka GFF files found in: {annotations_input}")

        gff_files = self._select_gff_files_for_binning(
            gff_files,
            binning_input,
            method_selection=method_selection,
        )

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

    def _select_gff_files_for_binning(
        self,
        gff_files: list[Path],
        binning_input: str | Path | list[str | Path] | tuple[str | Path, ...] | None,
        method_selection: set[str] | None = None,
    ) -> list[Path]:
        if not gff_files:
            return gff_files

        selected_methods = set(method_selection or set())
        if not selected_methods:
            selected_methods = detect_methods_from_inputs(binning_input)
        if not selected_methods:
            return gff_files

        methods_present: set[str] = set()
        filtered: list[Path] = []
        for gff_path in gff_files:
            method = _detect_method_from_text(str(gff_path))
            if method is not None:
                methods_present.add(method)
                if method in selected_methods:
                    filtered.append(gff_path)

        if not methods_present:
            return gff_files

        if filtered:
            logger.info(
                "Matched %d/%d Prokka GFF files to methods: %s.",
                len(filtered),
                len(gff_files),
                ", ".join(sorted(selected_methods)),
            )
            return filtered

        available_methods = ", ".join(sorted(methods_present))
        requested_methods = ", ".join(sorted(selected_methods))
        raise ValueError(
            "No Prokka annotations matched requested method(s) "
            f"'{requested_methods}' from {binning_input}. "
            f"Detected annotation methods in directory: {available_methods}. "
            "Provide annotations generated from the same binning method or pass "
            "a curated contig-level table via --annotation-data."
        )

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
    binning_input: str | Path | list[str | Path] | tuple[str | Path, ...] | None = None,
    method_selection: set[str] | None = None,
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
            binning_input=binning_input,
            method_selection=method_selection,
        )
    raise FileNotFoundError(f"Annotations input not found: {annotations_input}")
