import csv
import json
import re
import time
from pathlib import Path


SPIKE_START = 21563
SPIKE_END = 25384
UNKNOWN_AMINO_ACIDS = {"", "X", "_", "-"}

BASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = BASE_DIR.parent.parent
RAW_DATA_DIR = REPO_ROOT / "PRIEST" / "src" / "PRIEST_data" / "Raw Data"
OUTPUT_PATH = BASE_DIR / "priest_site_scores.csv"


def load_reference_spike():
    """Load and translate the reference Spike sequence.

    Args:
        None.

    Returns:
        The reference Spike amino-acid sequence without the terminal stop codon.
    """
    with (BASE_DIR / "genome.txt").open() as handle:
        next(handle)
        genome = "".join(line.strip() for line in handle)

    with (BASE_DIR / "codon_aa_mapping.json").open() as handle:
        codon_mapping = json.load(handle)

    spike_nt = genome[SPIKE_START - 1 : SPIKE_END]
    amino_acids = []
    for index in range(0, len(spike_nt) - 2, 3):
        amino_acids.append(codon_mapping.get(spike_nt[index : index + 3].upper(), "X"))

    reference = "".join(amino_acids)
    return reference[:-1] if reference.endswith("*") else reference


def iter_period_files():
    """Yield PRIEST raw-data files grouped by quarter.

    Args:
        None.

    Returns:
        An iterator of ``(period, [csv_paths])`` pairs for all discovered quarters.
    """
    quarter_raw_dir = RAW_DATA_DIR / "quarter_raw"
    if quarter_raw_dir.exists():
        for directory in sorted(quarter_raw_dir.iterdir()):
            if not directory.is_dir():
                continue

            match = re.fullmatch(r"year_(\d{4})_(\d)", directory.name)
            if not match:
                continue

            year = int(match.group(1))
            quarter = int(match.group(2)) + 1
            files = sorted(path for path in directory.glob("*.csv") if path.is_file())
            if files:
                yield f"Q{quarter}-{year}", files

    data_2023_dir = RAW_DATA_DIR / "2023 data"
    if data_2023_dir.exists():
        for file_path in sorted(data_2023_dir.glob("quarter_*_data.csv")):
            match = re.fullmatch(r"quarter_(\d+)_data\.csv", file_path.name)
            if match:
                yield f"Q{int(match.group(1))}-2023", [file_path]


def period_sort_key(period):
    """Create a chronological sort key for quarter labels.

    Args:
        period: Quarter label such as ``Q3-2021``.

    Returns:
        A tuple sortable by year and quarter.
    """
    quarter, year = period.split("-")
    return (int(year), int(quarter[1:]))


def build_lookup_rows():
    """Compute period-specific and global PRIEST site-score rows.

    Args:
        None.

    Returns:
        A list of dictionaries with ``period``, ``aa_position``, and
        ``priest_score`` ready to write to CSV.
    """
    reference_spike = load_reference_spike()
    global_mutated = [0] * (len(reference_spike) + 1)
    global_observed = [0] * (len(reference_spike) + 1)
    rows = []

    for period, files in sorted(iter_period_files(), key=lambda item: period_sort_key(item[0])):
        mutated = [0] * (len(reference_spike) + 1)
        observed = [0] * (len(reference_spike) + 1)
        sequence_count = 0

        for file_path in files:
            with file_path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    sequence = (row.get("Sequence") or "").strip().upper()
                    if not sequence:
                        continue

                    sequence_count += 1
                    limit = min(len(sequence), len(reference_spike))
                    for index in range(limit):
                        alt_aa = sequence[index]
                        if alt_aa in UNKNOWN_AMINO_ACIDS:
                            continue

                        position = index + 1
                        observed[position] += 1
                        global_observed[position] += 1

                        if alt_aa == reference_spike[index]:
                            continue

                        mutated[position] += 1
                        global_mutated[position] += 1

        for position in range(1, len(reference_spike) + 1):
            if observed[position]:
                rows.append(
                    {
                        "period": period,
                        "aa_position": position,
                        "priest_score": round(mutated[position] / observed[position], 6),
                    }
                )

        print(period, "sequences", sequence_count)

    for position in range(1, len(reference_spike) + 1):
        if global_observed[position]:
            rows.append(
                {
                    "period": "",
                    "aa_position": position,
                    "priest_score": round(global_mutated[position] / global_observed[position], 6),
                }
            )

    return rows


def main():
    """Build and write the local PRIEST lookup CSV.

    Args:
        None.

    Returns:
        None. Writes ``priest_site_scores.csv`` and prints summary stats.
    """
    start = time.time()
    rows = build_lookup_rows()
    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["period", "aa_position", "priest_score"])
        writer.writeheader()
        writer.writerows(rows)

    print("wrote", OUTPUT_PATH)
    print("rows", len(rows))
    print("seconds", round(time.time() - start, 2))


if __name__ == "__main__":
    main()
