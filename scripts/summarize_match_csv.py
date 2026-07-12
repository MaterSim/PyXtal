#!/usr/bin/env python3
"""Summarize ``yes`` counts in MLP-relax CSV outputs per case folder.

Counts matches in the ``original_matched`` and ``relaxed_matched`` columns,
and the number of unique relaxation attempts (``status`` not
``not_selected`` or ``propagated``), written by ``example_10_mlp_relax.py``.

Example:
    python summarize_match_csv.py Tests-0611
    python summarize_match_csv.py Tests-0611 --pattern 'MACEOFF_energy_15.0.csv'
    python summarize_match_csv.py Tests-0611 --sort attempted --min-total 10
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


RELAX_ATTEMPTED_STATUSES = frozenset({
    "ok",
    "relax_failed",
    "stall_skipped",
    "read_error",
    "relax_exception",
})


def count_csv_stats(csv_path: Path) -> tuple[int, int, int]:
    """Return (original_matched yes, relaxed_matched yes, relax attempts)."""
    orig = 0
    relax = 0
    attempted = 0
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return 0, 0, 0
        for row in reader:
            status = row.get("status", "").strip().lower()
            if status in RELAX_ATTEMPTED_STATUSES or (
                status and status not in {"not_selected", "propagated"}
            ):
                attempted += 1
            if row.get("original_matched", "").strip().lower() == "yes":
                orig += 1
            if row.get("relaxed_matched", "").strip().lower() == "yes":
                relax += 1
    return orig, relax, attempted


def collect_case_stats(root: Path, pattern: str) -> list[dict]:
    """Scan case subfolders and aggregate yes counts per case."""
    rows = []
    for case_dir in sorted(root.iterdir()):
        if not case_dir.is_dir():
            continue

        case_orig = 0
        case_relax = 0
        case_attempted = 0
        matched_files: list[str] = []
        found_csv = False

        for csv_path in sorted(case_dir.glob(pattern)):
            if not csv_path.is_file():
                continue
            found_csv = True
            orig, relax, attempted = count_csv_stats(csv_path)
            case_orig += orig
            case_relax += relax
            case_attempted += attempted
            if orig or relax:
                matched_files.append(f"{csv_path.name}:{orig}+{relax}")

        if not found_csv:
            continue

        rows.append(
            {
                "case": case_dir.name,
                "orig": case_orig,
                "relax": case_relax,
                "attempted": case_attempted,
                "total": case_orig + case_relax,
                "files": ", ".join(matched_files),
            }
        )
    return rows


def print_table(rows: list[dict], sort_key: str) -> None:
    header = f"{'CASE':<14} {'orig':>6} {'relax':>6} {'attempt':>7}"
    print(header)
    print("-" * len(header))

    key_map = {
        "case": lambda r: r["case"],
        "orig": lambda r: (r["orig"], r["case"]),
        "relax": lambda r: (r["relax"], r["case"]),
        "attempted": lambda r: (r["attempted"], r["case"]),
        "total": lambda r: (r["total"], r["case"]),
    }
    for row in sorted(rows, key=key_map[sort_key], reverse=sort_key != "case"):
        print(
            f"{row['case']:<14} {row['orig']:>6} {row['relax']:>6} "
            f"{row['attempted']:>7}"
        )

    n_cases = len(rows)
    n_with_matches = sum(1 for r in rows if r["total"] > 0)
    sum_orig = sum(r["orig"] for r in rows)
    sum_relax = sum(r["relax"] for r in rows)
    sum_attempted = sum(r["attempted"] for r in rows)
    print("-" * len(header))
    print(
        f"{'TOTAL':<14} {sum_orig:>6} {sum_relax:>6} {sum_attempted:>7}  "
        f"({n_with_matches}/{n_cases} cases with matches)"
    )


def write_csv(rows: list[dict], out_path: Path, sort_key: str) -> None:
    key_map = {
        "case": lambda r: r["case"],
        "orig": lambda r: (r["orig"], r["case"]),
        "relax": lambda r: (r["relax"], r["case"]),
        "attempted": lambda r: (r["attempted"], r["case"]),
        "total": lambda r: (r["total"], r["case"]),
    }
    sorted_rows = sorted(rows, key=key_map[sort_key], reverse=sort_key != "case")
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["case", "orig", "relax", "attempted", "files"]
        )
        writer.writeheader()
        for row in sorted_rows:
            writer.writerow(
                {
                    "case": row["case"],
                    "orig": row["orig"],
                    "relax": row["relax"],
                    "attempted": row["attempted"],
                    "files": row["files"],
                }
            )
    print(f"\nWrote summary CSV: {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Count yes values in original_matched / relaxed_matched columns.",
    )
    parser.add_argument(
        "root",
        nargs="?",
        default="Tests-0611",
        help="Root directory containing per-case subfolders (default: Tests-0611)",
    )
    parser.add_argument(
        "--pattern",
        default="MACEOFF_energy_15.0.csv",
        help="CSV filename glob inside each case folder (default: MACEOFF_energy_15.0.csv)",
    )
    parser.add_argument(
        "--sort",
        choices=("case", "orig", "relax", "attempted", "total"),
        default="relax",
        help="Sort key (default: relax, descending)",
    )
    parser.add_argument(
        "--min-total",
        type=int,
        default=0,
        help="Only show cases with total yes count >= this value",
    )
    parser.add_argument(
        "--out-csv",
        default=None,
        help="Optional path to write the summary table as CSV",
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    if not root.is_dir():
        raise SystemExit(f"Not a directory: {root}")

    rows = collect_case_stats(root, args.pattern)
    if args.min_total > 0:
        rows = [r for r in rows if r["total"] >= args.min_total]

    print(f"Root: {root}")
    print(f"Pattern: {args.pattern}")
    print_table(rows, args.sort)

    if args.out_csv:
        write_csv(rows, Path(args.out_csv).resolve(), args.sort)


if __name__ == "__main__":
    main()
