"""Organism registry for the generic virus mutation platform.

Two paths are supported:

* **Built-in organisms** live under ``<project_root>/organisms/<name>/``.
  COVID and Influenza ship with the platform. Loaders cache the parsed result
  per-process so repeated predictions don't re-read the FASTA.

* **Custom organisms** live inside the user's bundle as helper files. The
  bundle's ``bundle_metadata.json`` declares ``genome_file`` (and optionally
  ``protein_regions_file``); the runtime resolves those paths inside
  ``uploaded_models/<bundle>/`` and reads them here.

Protein-region CSV format (``name,start,end``) uses **1-based inclusive**
coordinates — the same convention COVID's existing ``PROTEIN_REGIONS`` dict
already uses, so the rest of the pipeline keeps working without translation.
"""

import csv
import json
import os
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

GENOME_EXTRACTOR_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ORGANISMS_DIR = os.path.join(GENOME_EXTRACTOR_DIR, "organisms")

ProteinRegions = Dict[str, List[int]]


def _read_fasta_body(path: str) -> str:
    """Return the concatenated sequence from a FASTA file (headers skipped).

    Accepts multi-line FASTA and ignores blank lines. Matches the existing
    ``helpers.read_genome_sequence`` behavior so swapping between them is
    transparent.
    """
    parts = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith(">"):
                continue
            parts.append(stripped)
    return "".join(parts)


def _read_protein_regions_csv(path: str) -> ProteinRegions:
    """Parse a ``name,start,end`` CSV (1-based inclusive) into a dict.

    Returns ``{}`` when the file is missing — protein regions are optional.
    """
    regions: ProteinRegions = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            name = (row.get("name") or "").strip()
            if not name:
                continue
            try:
                start = int(row["start"])
                end = int(row["end"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid row in {path}: {row!r} — expected name,start,end"
                ) from exc
            if end <= start:
                raise ValueError(
                    f"Invalid protein region {name!r} in {path}: end must be > start"
                )
            regions[name] = [start, end]
    return regions


@lru_cache(maxsize=8)
def read_builtin_organism(name: str) -> Tuple[str, ProteinRegions]:
    """Load a built-in organism by name. Returns ``(genome_seq, protein_regions)``.

    Raises ``FileNotFoundError`` if the organism is not registered or the
    genome.fasta is missing. ``protein_regions`` may be empty when the organism
    has no annotated proteins.
    """
    organism_dir = os.path.join(ORGANISMS_DIR, name)
    if not os.path.isdir(organism_dir):
        raise FileNotFoundError(
            f"Built-in organism not found: {name!r}. "
            f"Available: {list_builtin_organisms()}"
        )

    genome_path = os.path.join(organism_dir, "genome.fasta")
    if not os.path.isfile(genome_path):
        raise FileNotFoundError(f"Missing genome.fasta in {organism_dir}")

    genome_seq = _read_fasta_body(genome_path)

    regions_path = os.path.join(organism_dir, "protein_regions.csv")
    protein_regions: ProteinRegions = {}
    if os.path.isfile(regions_path):
        protein_regions = _read_protein_regions_csv(regions_path)

    return genome_seq, protein_regions


def read_custom_organism_from_bundle(
    bundle_dir: str,
    genome_file: str,
    protein_regions_file: Optional[str] = None,
) -> Tuple[str, ProteinRegions]:
    """Load a custom organism shipped inside a bundle's helper files.

    Paths are resolved relative to ``bundle_dir``; the helper file names must
    have already passed sanitization. Raises ``FileNotFoundError`` if the
    declared genome file is not present.
    """
    if not genome_file:
        raise ValueError("custom organism bundle must declare genome_file")

    genome_path = os.path.join(bundle_dir, genome_file)
    if not os.path.isfile(genome_path):
        raise FileNotFoundError(
            f"Custom genome file not found in bundle: {genome_file} "
            f"(looked in {bundle_dir})"
        )

    genome_seq = _read_fasta_body(genome_path)

    protein_regions: ProteinRegions = {}
    if protein_regions_file:
        regions_path = os.path.join(bundle_dir, protein_regions_file)
        if os.path.isfile(regions_path):
            protein_regions = _read_protein_regions_csv(regions_path)

    return genome_seq, protein_regions


def read_organism_metadata(name: str) -> Dict:
    """Load ``organism.json`` for a built-in organism, or ``{}`` if absent.

    Used for surfacing display names / accessions to the UI; not required for
    the predict pipeline itself.
    """
    organism_dir = os.path.join(ORGANISMS_DIR, name)
    metadata_path = os.path.join(organism_dir, "organism.json")
    if not os.path.isfile(metadata_path):
        return {}
    with open(metadata_path, "r", encoding="utf-8") as handle:
        loaded = json.load(handle)
    return loaded if isinstance(loaded, dict) else {}


def list_builtin_organisms() -> List[str]:
    """Names of all built-in organisms (alphabetical)."""
    if not os.path.isdir(ORGANISMS_DIR):
        return []
    return sorted(
        name
        for name in os.listdir(ORGANISMS_DIR)
        if os.path.isdir(os.path.join(ORGANISMS_DIR, name)) and not name.startswith(".")
    )


# ---------------------------------------------------------------------------
# Variant catalog
# ---------------------------------------------------------------------------
# Built-in organisms can ship a curated catalog of variant FASTAs under
# organisms/<name>/variants/, indexed by variants.csv. The reference sequence
# (genome.fasta) and protein_regions.csv stay authoritative; variants are
# strict overrides for the genome only. This matches the biological reality
# for influenza HA where each variant is a distinct sequence (not a set of
# point mutations) and where the reference's protein_regions are the right
# annotation surface to overlay against.


@lru_cache(maxsize=16)
def list_variants(organism_name: str) -> List[Dict[str, str]]:
    """Return the catalog of variants for a built-in organism, or ``[]``.

    Each entry is a dict with the columns declared in ``variants.csv``:
    ``name``, ``display_name``, ``fasta_file``, ``accession``, ``cds_start``,
    ``cds_end``, ``year``, ``subtype``, ``description``. The frontend uses
    these to drive the variant dropdown.
    """
    variants_path = os.path.join(ORGANISMS_DIR, organism_name, "variants", "variants.csv")
    if not os.path.isfile(variants_path):
        return []

    catalog: List[Dict[str, str]] = []
    with open(variants_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            name = (row.get("name") or "").strip()
            fasta_file = (row.get("fasta_file") or "").strip()
            if not name or not fasta_file:
                continue
            catalog.append({key: (value or "").strip() for key, value in row.items()})
    return catalog


INFLUENZA_SUBTYPE_ORGANISMS = ("influenza_h1n1", "influenza_h3n2", "influenza_h5n1")


@lru_cache(maxsize=1)
def list_all_influenza_variants() -> List[Dict[str, str]]:
    """Return every influenza HA strain across all three built-in subtypes as
    a flat catalog.

    Each entry carries an extra ``organism`` field (e.g. ``influenza_h1n1``)
    pointing at the subtype folder that owns its FASTA and protein_regions.
    The frontend renders the entries in one dropdown; the backend uses
    ``find_variant_subtype()`` to dispatch a generic ``organism="influenza"``
    bundle to the right subtype-specific protein_regions.

    Ordering: subtype groups (H1N1 → H3N2 → H5N1), references first within
    each subtype, then chronologically.
    """
    consolidated: List[Dict[str, str]] = []
    for organism in INFLUENZA_SUBTYPE_ORGANISMS:
        entries = list_variants(organism)
        # Stable ordering: reference rows first, then by year ascending.
        entries = sorted(
            entries,
            key=lambda row: (
                0 if str(row.get("is_reference", "")).lower() == "true" else 1,
                int(row.get("year") or 0),
            ),
        )
        for entry in entries:
            consolidated.append({**entry, "organism": organism})
    return consolidated


@lru_cache(maxsize=32)
def find_variant_subtype(variant_name: str) -> Optional[str]:
    """Return the subtype organism that owns ``variant_name``, or ``None``.

    Used when a bundle declares the generic ``organism="influenza"`` and the
    runtime needs to pick the right subtype-specific ``protein_regions.csv``.
    """
    for entry in list_all_influenza_variants():
        if entry["name"] == variant_name:
            return entry["organism"]
    return None


@lru_cache(maxsize=32)
def read_variant(organism_name: str, variant_name: str) -> str:
    """Load a variant's genome sequence as a string.

    Raises ``FileNotFoundError`` if the variant is not in the organism's
    catalog or the referenced FASTA file is missing.
    """
    catalog = list_variants(organism_name)
    entry = next((item for item in catalog if item["name"] == variant_name), None)
    if entry is None:
        available = [item["name"] for item in catalog]
        raise FileNotFoundError(
            f"Variant {variant_name!r} not in {organism_name}'s catalog. "
            f"Available: {available}"
        )

    variant_path = os.path.join(
        ORGANISMS_DIR, organism_name, "variants", entry["fasta_file"]
    )
    if not os.path.isfile(variant_path):
        raise FileNotFoundError(
            f"Variant FASTA missing on disk: {variant_path}"
        )
    return _read_fasta_body(variant_path)


def get_variant_mutations(organism: str, node_id) -> list:
    """Return the parsed mutation list for organism-specific feature extraction.

    Only COVID uses a node-id-based mutation diff table; all other organisms
    return an empty list.
    """
    if organism == "covid":
        from .feature_extractor_updated import parse_mutations
        return parse_mutations(node_id) if node_id else []
    return []


def apply_variant_sequence(organism: str, genome_seq: str, node_id) -> str:
    """Apply organism-specific variant overlay and return the modified genome.

    Only COVID uses a node-id-based mutation diff table. All other organisms
    return the genome sequence unchanged.
    """
    if organism == "covid":
        from .feature_extractor_updated import construct_variant_genome
        mutations = get_variant_mutations(organism, node_id)
        return construct_variant_genome(genome_seq, mutations)
    return genome_seq
