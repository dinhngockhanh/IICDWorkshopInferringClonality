#!/usr/bin/env python3
"""Build one canonical S8 table, then the five algorithm input files.

Raw inputs (read only):

- data/raw/S8-noXY_mutation_table.csv
- data/raw/S8-noXY.battenberg.txt
- data/raw/S8-noXY.cellularity_ploidy.txt

Battenberg solution A is the only copy-number solution used. It matches
major_cn/minor_cn on every mutation that falls inside a segment. Mutations
outside those segments keep the copy number already stored on the mutation.
Non-diploid ones are wrapped in extra segments so PhyloWGS does not treat
them as diploid. Homozygous deletions (0, 0) are kept in the canonical table
but left out of PyClone-VI and out of PhyloWGS CNV records.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PREPARED = ROOT / "data" / "prepared"
SAMPLE_ID = "S8-noXY"
MU_R = 0.999
MU_V = 0.5


def normalize_chrom(value: str) -> str:
    text = str(value).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    if text.upper() == "MT":
        return "M"
    return text


def chrom_sort_key(chrom: str):
    if chrom.isdigit():
        return (0, int(chrom), chrom)
    named = {"X": 23, "Y": 24, "M": 25}
    if chrom.upper() in named:
        return (1, named[chrom.upper()], chrom)
    return (2, 0, chrom)


def read_sample(path: Path) -> tuple[float, float]:
    lines = path.read_text().splitlines()
    if not lines:
        sys.exit(f"Empty purity file: {path}")
    delimiter = "\t" if "\t" in lines[0] else "," if "," in lines[0] else None
    if delimiter is None:
        header = lines[0].split()
        values = lines[1].split()
        row = dict(zip(header, values))
    else:
        row = next(csv.DictReader(lines, delimiter=delimiter))
    return float(row["cellularity"]), float(row["ploidy"])


def read_mutations(path: Path) -> list[dict]:
    mutations = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            chrom = normalize_chrom(row["Chromosome"])
            pos = int(float(row["Position"]))
            ref = int(row["Ref_count"])
            alt = int(row["Alt_count"])
            major = int(row["major_cn"])
            minor = int(row["minor_cn"])
            multiplicity = int(row["multiplicity"])
            depth = ref + alt
            mutations.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "ref_count": ref,
                    "alt_count": alt,
                    "major_cn": major,
                    "minor_cn": minor,
                    "multiplicity": multiplicity,
                    "depth": depth,
                    "vaf": (alt / depth) if depth else None,
                }
            )
    mutations.sort(key=lambda item: (chrom_sort_key(item["chrom"]), item["pos"]))
    for item in mutations:
        item["mutation_id"] = f"{item['chrom']}:{item['pos']}"
    ids = [item["mutation_id"] for item in mutations]
    if len(ids) != len(set(ids)):
        sys.exit("Mutation ids are not unique; chrom:pos collided.")
    return mutations


def _parse_cn(value: str):
    if value is None or value == "" or value.upper() == "NA":
        return None
    return int(float(value))


def read_battenberg(path: Path) -> tuple[list[dict], list[str]]:
    segments = []
    skipped = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            major = _parse_cn(row["nMaj1_A"])
            minor = _parse_cn(row["nMin1_A"])
            if major is None or minor is None:
                skipped.append(f"{row['chr']}:{row['startpos']}-{row['endpos']}")
                continue
            frac_text = row["frac1_A"]
            if frac_text in ("", "NA"):
                frac = 1.0
            else:
                frac = float(frac_text)
            segments.append(
                {
                    "chrom": normalize_chrom(row["chr"]),
                    "start": int(float(row["startpos"])),
                    "end": int(float(row["endpos"])),
                    "major_cn": major,
                    "minor_cn": minor,
                    "frac": frac,
                    "source": "battenberg",
                }
            )
    segments.sort(key=lambda item: (chrom_sort_key(item["chrom"]), item["start"], item["end"]))
    return segments, skipped


def segment_contains(segment: dict, mutation: dict) -> bool:
    return (
        segment["chrom"] == mutation["chrom"]
        and segment["start"] <= mutation["pos"] <= segment["end"]
    )


def interval_is_blocked(chrom: str, left: int, right: int, segments: list[dict], mutations: list[dict]) -> bool:
    for segment in segments:
        if segment["chrom"] != chrom:
            continue
        if segment["start"] < right and segment["end"] > left:
            return True
    for mutation in mutations:
        if mutation["chrom"] != chrom:
            continue
        if left < mutation["pos"] < right:
            return True
    return False


def gap_fill_segments(mutations: list[dict], battenberg: list[dict]) -> list[dict]:
    """Wrap uncovered non-diploid mutations so PhyloWGS sees their copy number.

    Homozygous deletions and diploid mutations are not wrapped. A run is merged
    only when the open interval between its mutations contains neither another
    mutation nor a Battenberg segment.
    """
    covered = {
        mutation["mutation_id"]
        for mutation in mutations
        if any(segment_contains(segment, mutation) for segment in battenberg)
    }
    candidates = [
        mutation
        for mutation in mutations
        if mutation["mutation_id"] not in covered
        and not (mutation["major_cn"] == 0 and mutation["minor_cn"] == 0)
        and not (mutation["major_cn"] == 1 and mutation["minor_cn"] == 1)
    ]
    groups: list[list[dict]] = []
    for mutation in candidates:
        if not groups:
            groups.append([mutation])
            continue
        previous = groups[-1][-1]
        same_state = (
            mutation["chrom"] == previous["chrom"]
            and mutation["major_cn"] == previous["major_cn"]
            and mutation["minor_cn"] == previous["minor_cn"]
        )
        blocked = interval_is_blocked(
            mutation["chrom"], previous["pos"], mutation["pos"], battenberg, mutations
        )
        if same_state and not blocked:
            groups[-1].append(mutation)
        else:
            groups.append([mutation])

    filled = []
    for group in groups:
        filled.append(
            {
                "chrom": group[0]["chrom"],
                "start": group[0]["pos"],
                "end": group[-1]["pos"],
                "major_cn": group[0]["major_cn"],
                "minor_cn": group[0]["minor_cn"],
                "frac": 1.0,
                "source": "gap_fill",
            }
        )
    return filled


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_csv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def is_diploid(mutation: dict) -> bool:
    return mutation["major_cn"] == 1 and mutation["minor_cn"] == 1


def is_zero_cn(mutation: dict) -> bool:
    return mutation["major_cn"] == 0 and mutation["minor_cn"] == 0


def total_cn(record: dict) -> int:
    return record["major_cn"] + record["minor_cn"]


def check_copy_number_coverage(mutations: list[dict], segments: list[dict]) -> None:
    for mutation in mutations:
        covering = [segment for segment in segments if segment_contains(segment, mutation)]
        if is_diploid(mutation):
            bad = [segment for segment in covering if total_cn(segment) != 2]
            if bad:
                sys.exit(f"Diploid mutation {mutation['mutation_id']} falls in a non-diploid segment.")
            continue
        if is_zero_cn(mutation):
            bad = [segment for segment in covering if total_cn(segment) != 0]
            if bad:
                sys.exit(f"Zero-CN mutation {mutation['mutation_id']} falls in a non-zero segment.")
            continue
        matched = [
            segment
            for segment in covering
            if segment["major_cn"] == mutation["major_cn"] and segment["minor_cn"] == mutation["minor_cn"]
        ]
        if len(matched) != 1:
            sys.exit(
                f"Mutation {mutation['mutation_id']} CN "
                f"({mutation['major_cn']},{mutation['minor_cn']}) is covered by {len(matched)} matching segments."
            )


def main() -> None:
    mutation_path = RAW / "S8-noXY_mutation_table.csv"
    battenberg_path = RAW / "S8-noXY.battenberg.txt"
    purity_path = RAW / "S8-noXY.cellularity_ploidy.txt"
    for path in (mutation_path, battenberg_path, purity_path):
        if not path.exists():
            sys.exit(f"Missing raw input: {path}")

    purity, ploidy = read_sample(purity_path)
    mutations = [item for item in read_mutations(mutation_path) if item["depth"] > 0]
    battenberg, skipped_segments = read_battenberg(battenberg_path)
    gap_fill = gap_fill_segments(mutations, battenberg)
    segments = battenberg + gap_fill
    segments.sort(key=lambda item: (chrom_sort_key(item["chrom"]), item["start"], item["end"], item["source"]))
    check_copy_number_coverage(mutations, segments)

    write_tsv(
        PREPARED / "mutations.tsv",
        [
            "mutation_id",
            "chrom",
            "pos",
            "ref_count",
            "alt_count",
            "major_cn",
            "minor_cn",
            "multiplicity",
            "depth",
            "vaf",
        ],
        [
            [
                item["mutation_id"],
                item["chrom"],
                item["pos"],
                item["ref_count"],
                item["alt_count"],
                item["major_cn"],
                item["minor_cn"],
                item["multiplicity"],
                item["depth"],
                f"{item['vaf']:.12g}",
            ]
            for item in mutations
        ],
    )
    write_tsv(
        PREPARED / "segments.tsv",
        ["chrom", "start", "end", "major_cn", "minor_cn", "frac", "source"],
        [
            [
                item["chrom"],
                item["start"],
                item["end"],
                item["major_cn"],
                item["minor_cn"],
                f"{item['frac']:.10g}",
                item["source"],
            ]
            for item in segments
        ],
    )
    write_tsv(
        PREPARED / "sample.tsv",
        ["purity", "ploidy"],
        [[f"{purity:.16g}", f"{ploidy:.16g}"]],
    )

    diploid = [
        item
        for item in mutations
        if is_diploid(item) and item["multiplicity"] == 1 and item["vaf"] is not None
    ]
    decode_rows = [item for item in mutations if is_diploid(item)]
    mobster_rows = [
        item for item in diploid if item["vaf"] is not None and 0 < item["vaf"] < 1
    ]
    pyclone_rows = [
        item for item in mutations if item["multiplicity"] == 1 and total_cn(item) != 0
    ]
    zero_cn_ids = [item["mutation_id"] for item in mutations if is_zero_cn(item)]

    inputs = PREPARED / "inputs"
    write_tsv(
        inputs / "decode" / "mutations.tsv",
        ["mutation_id", "chrom", "pos", "Ref_count", "Alt_count"],
        [
            [item["mutation_id"], item["chrom"], item["pos"], item["ref_count"], item["alt_count"]]
            for item in decode_rows
        ],
    )
    write_csv(
        inputs / "mobster" / f"{SAMPLE_ID}.vaf.csv",
        ["mutation_id", "vaf", "alt_count", "depth"],
        [
            [
                item["mutation_id"],
                f"{item['alt_count'] / item['depth']:.6f}",
                item["alt_count"],
                item["depth"],
            ]
            for item in mobster_rows
        ],
    )

    sciclone_dir = inputs / "sciclone"
    write_tsv(
        sciclone_dir / "vafs.tsv",
        ["chr", "pos", "ref_reads", "var_reads", "vaf"],
        [
            [
                item["chrom"],
                item["pos"],
                item["ref_count"],
                item["alt_count"],
                f"{100.0 * item['alt_count'] / item['depth']:.4f}",
            ]
            for item in mobster_rows
        ],
    )
    write_tsv(
        sciclone_dir / "cn.tsv",
        ["chr", "start", "stop", "segment_mean"],
        [
            [item["chrom"], item["start"], item["end"], total_cn(item)]
            for item in segments
        ],
    )
    loh_path = sciclone_dir / "exclude.loh.bed"
    loh_path.parent.mkdir(parents=True, exist_ok=True)
    with loh_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for item in segments:
            if item["minor_cn"] == 0:
                writer.writerow([item["chrom"], item["start"], item["end"]])

    write_tsv(
        inputs / "pyclonevi" / f"{SAMPLE_ID}.tsv",
        [
            "mutation_id",
            "sample_id",
            "ref_counts",
            "alt_counts",
            "major_cn",
            "minor_cn",
            "normal_cn",
            "tumour_content",
        ],
        [
            [
                item["mutation_id"],
                SAMPLE_ID,
                item["ref_count"],
                item["alt_count"],
                item["major_cn"],
                item["minor_cn"],
                2,
                f"{purity:.16g}",
            ]
            for item in pyclone_rows
        ],
    )

    for index, item in enumerate(mutations):
        item["sid"] = f"s{index}"
    phylowgs_dir = inputs / "phylowgs"
    write_tsv(
        phylowgs_dir / "ssm_data.txt",
        ["id", "gene", "a", "d", "mu_r", "mu_v"],
        [
            [item["sid"], f"chr{item['chrom']}:{item['pos']}", item["ref_count"], item["depth"], MU_R, MU_V]
            for item in mutations
        ],
    )
    write_tsv(
        phylowgs_dir / "ssm_map.tsv",
        ["sid", "mutation_id", "chrom", "pos"],
        [[item["sid"], item["mutation_id"], item["chrom"], item["pos"]] for item in mutations],
    )

    cnv_segments = []
    for segment in segments:
        if total_cn(segment) == 2 or total_cn(segment) == 0:
            continue
        inside = [item for item in mutations if segment_contains(segment, item)]
        if not inside:
            continue
        cnv_segments.append((segment, inside))
    with (phylowgs_dir / "cnv_data.txt").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["cnv", "a", "d", "ssms", "physical_cnvs"])
        for index, (segment, inside) in enumerate(cnv_segments):
            ref_sum = sum(item["ref_count"] for item in inside)
            depth_sum = sum(item["depth"] for item in inside)
            ssms = ";".join(f"{item['sid']},{item['multiplicity']},0" for item in inside)
            physical = (
                f"chrom={segment['chrom']},start={segment['start']},end={segment['end']},"
                f"major_cn={segment['major_cn']},minor_cn={segment['minor_cn']},"
                f"cell_prev={segment['frac']:.10g}"
            )
            writer.writerow([f"c{index}", ref_sum, depth_sum, ssms, physical])

    uncovered = [
        item
        for item in mutations
        if not any(segment_contains(segment, item) for segment in battenberg)
    ]
    cn_counts = Counter(f"{item['major_cn']},{item['minor_cn']}" for item in mutations)
    mult_counts = Counter(str(item["multiplicity"]) for item in mutations)
    subclonal = [item for item in battenberg if abs(item["frac"] - 1.0) > 1e-8]
    qc = {
        "sample_id": SAMPLE_ID,
        "purity": purity,
        "ploidy": ploidy,
        "n_mutations": len(mutations),
        "copy_number_counts": dict(cn_counts),
        "multiplicity_counts": dict(mult_counts),
        "n_battenberg_segments": len(battenberg),
        "n_battenberg_segments_skipped": len(skipped_segments),
        "skipped_battenberg_segments": skipped_segments,
        "n_subclonal_battenberg_segments": len(subclonal),
        "subclonal_battenberg_segments": [
            {
                "chrom": item["chrom"],
                "start": item["start"],
                "end": item["end"],
                "major_cn": item["major_cn"],
                "minor_cn": item["minor_cn"],
                "frac": item["frac"],
            }
            for item in subclonal
        ],
        "n_uncovered_mutations": len(uncovered),
        "n_gap_fill_mutations": 0,
        "n_gap_fill_segments": len(gap_fill),
        "gap_fill_segments": [
            {
                "chrom": item["chrom"],
                "start": item["start"],
                "end": item["end"],
                "major_cn": item["major_cn"],
                "minor_cn": item["minor_cn"],
            }
            for item in gap_fill
        ],
        "dropped_zero_cn_mutation_ids": zero_cn_ids,
        "n_decode": len(decode_rows),
        "n_mobster": len(mobster_rows),
        "n_sciclone": len(mobster_rows),
        "n_pyclonevi": len(pyclone_rows),
        "n_phylowgs_ssm": len(mutations),
        "n_phylowgs_cnv": len(cnv_segments),
        "n_sciclone_depth_below_10": sum(1 for item in mobster_rows if item["depth"] < 10),
    }
    # Recount gap-fill mutations directly. The expression above is hard to read.
    gap_ids = []
    for segment in gap_fill:
        for item in mutations:
            if segment_contains(segment, item) and not is_diploid(item) and not is_zero_cn(item):
                gap_ids.append(item["mutation_id"])
    qc["n_gap_fill_mutations"] = len(set(gap_ids))
    qc["gap_fill_mutation_ids"] = sorted(set(gap_ids), key=lambda value: (chrom_sort_key(value.split(":")[0]), int(value.split(":")[1])))

    expected = {
        "n_mutations": 10751,
        "n_battenberg_segments": 27,
        "n_subclonal_battenberg_segments": 3,
        "n_uncovered_mutations": 137,
        "n_gap_fill_mutations": 25,
        "n_decode": 7898,
        "n_mobster": 7898,
        "n_pyclonevi": 10496,
        "n_phylowgs_ssm": 10751,
        "n_zero_cn": 5,
    }
    actual = {
        "n_mutations": qc["n_mutations"],
        "n_battenberg_segments": qc["n_battenberg_segments"],
        "n_subclonal_battenberg_segments": qc["n_subclonal_battenberg_segments"],
        "n_uncovered_mutations": qc["n_uncovered_mutations"],
        "n_gap_fill_mutations": qc["n_gap_fill_mutations"],
        "n_decode": qc["n_decode"],
        "n_mobster": qc["n_mobster"],
        "n_pyclonevi": qc["n_pyclonevi"],
        "n_phylowgs_ssm": qc["n_phylowgs_ssm"],
        "n_zero_cn": len(zero_cn_ids),
    }
    mismatches = {key: (actual[key], expected[key]) for key in expected if actual[key] != expected[key]}
    if mismatches:
        sys.exit(f"Prepared counts do not match the S8 audit: {mismatches}")
    if cn_counts["1,1"] != 7898 or cn_counts["0,0"] != 5:
        sys.exit(f"Unexpected copy-number counts: {cn_counts}")

    qc_path = PREPARED / "qc.json"
    qc_path.write_text(json.dumps(qc, indent=2) + "\n")
    print(f"Wrote {PREPARED}")
    print(
        f"mutations={qc['n_mutations']} decode={qc['n_decode']} "
        f"mobster={qc['n_mobster']} pyclonevi={qc['n_pyclonevi']} "
        f"phylowgs_ssm={qc['n_phylowgs_ssm']} phylowgs_cnv={qc['n_phylowgs_cnv']} "
        f"gap_fill_segments={qc['n_gap_fill_segments']}"
    )


if __name__ == "__main__":
    main()
