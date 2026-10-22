#!/usr/bin/env python3
"""Join each tool's output onto the full S8 mutation table.

One assignment file is written per tool:

results/assignments/<sample>_<tool>.tsv

Columns are mutation_id, chrom, pos, tool, included, cluster, ccf, note.
A mutation the tool did not use has included=false.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ("decode", "mobster", "sciclone", "pyclonevi", "phylowgs")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sample_id(root: Path) -> str:
    sample = read_tsv(root / "data" / "prepared" / "sample.tsv")
    mutations = read_tsv(root / "data" / "prepared" / "mutations.tsv")
    if not mutations:
        sys.exit("data/prepared/mutations.tsv is empty. Run src/00_prepare.py first.")
    del sample
    config = root / "config" / "params.toml"
    for line in config.read_text().splitlines():
        if line.startswith("sample_id"):
            return line.split("=", 1)[1].strip().strip('"').strip("'")
    sys.exit("sample_id is missing from config/params.toml")


def write_assignments(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ["mutation_id", "chrom", "pos", "tool", "included", "cluster", "ccf", "note"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def base_rows(mutations: list[dict[str, str]], tool: str) -> list[dict[str, str]]:
    return [
        {
            "mutation_id": item["mutation_id"],
            "chrom": item["chrom"],
            "pos": item["pos"],
            "tool": tool,
            "included": "false",
            "cluster": "",
            "ccf": "",
            "note": "",
        }
        for item in mutations
    ]


def fill_from_components(
    rows: list[dict[str, str]],
    components: list[dict[str, str]],
    missing_note: str,
    input_ids: set[str] | None = None,
    filtered_note: str = "filtered after input",
) -> None:
    by_id = {row["mutation_id"]: row for row in components}
    if len(by_id) != len(components):
        sys.exit("Component table has duplicate mutation ids.")
    for row in rows:
        found = by_id.get(row["mutation_id"])
        if found is None:
            if input_ids is not None and row["mutation_id"] in input_ids:
                row["note"] = filtered_note
            else:
                row["note"] = missing_note
            continue
        row["included"] = "true"
        row["cluster"] = found.get("cluster", "")
        row["ccf"] = found.get("ccf", "")
        row["note"] = found.get("note", "")


def run_r(root: Path, tool: str) -> None:
    script = root / "src" / "extract_r.R"
    completed = subprocess.run(["Rscript", str(script), tool], cwd=root)
    if completed.returncode != 0:
        sys.exit(f"R extraction failed for {tool}")


def ids_from_column(path: Path, column: str, delimiter: str = "\t") -> set[str]:
    with path.open(newline="") as handle:
        return {row[column] for row in csv.DictReader(handle, delimiter=delimiter)}


def extract_r_tool(root: Path, tool: str, mutations: list[dict[str, str]], sample: str) -> None:
    run_r(root, tool)
    components = read_tsv(root / "results" / tool / "components.tsv")
    rows = base_rows(mutations, tool)
    if tool == "decode":
        fill_from_components(rows, components, "unassigned")
    elif tool == "mobster":
        input_ids = ids_from_column(
            root / "data" / "prepared" / "inputs" / "mobster" / f"{sample}.vaf.csv",
            "mutation_id",
            ",",
        )
        fill_from_components(rows, components, "not in tool input", input_ids, "missing from MOBSTER fit")
    elif tool == "sciclone":
        vafs = read_tsv(root / "data" / "prepared" / "inputs" / "sciclone" / "vafs.tsv")
        input_ids = {f"{item['chr']}:{item['pos']}" for item in vafs}
        fill_from_components(
            rows,
            components,
            "not in tool input",
            input_ids,
            "filtered by sciClone",
        )
    else:
        sys.exit(f"Unexpected R tool {tool}")
    write_assignments(root / "results" / "assignments" / f"{sample}_{tool}.tsv", rows)


def extract_pyclone(root: Path, mutations: list[dict[str, str]], sample: str) -> None:
    result_path = root / "results" / "pyclonevi" / f"{sample}_results.tsv"
    if not result_path.exists():
        sys.exit(f"Missing PyClone-VI results: {result_path}")
    results = read_tsv(result_path)
    components = [
        {
            "mutation_id": item["mutation_id"],
            "cluster": item["cluster_id"],
            "ccf": item["cellular_prevalence"],
            "note": "",
        }
        for item in results
    ]
    input_ids = ids_from_column(
        root / "data" / "prepared" / "inputs" / "pyclonevi" / f"{sample}.tsv",
        "mutation_id",
    )
    rows = base_rows(mutations, "pyclonevi")
    fill_from_components(rows, components, "not in tool input", input_ids, "missing from PyClone-VI results")
    write_assignments(root / "results" / "assignments" / f"{sample}_pyclonevi.tsv", rows)


def choose_tree(summary: dict) -> str:
    densities = summary.get("tree_densities") or {}
    if densities:
        return str(max(densities, key=lambda key: float(densities[key])))
    trees = summary.get("trees") or {}
    if not trees:
        sys.exit("PhyloWGS summary contains no trees.")
    return str(next(iter(trees)))


def load_summary(path: Path) -> dict:
    with gzip.open(path, "rt") as handle:
        return json.load(handle)


def load_mutass(path: Path, tree_id: str) -> dict:
    if not zipfile.is_zipfile(path):
        sys.exit(f"PhyloWGS mutass is not a zip file: {path}")
    candidates = [f"{tree_id}.json"]
    try:
        candidates.append(f"{int(float(tree_id))}.json")
    except ValueError:
        pass
    with zipfile.ZipFile(path) as archive:
        names = [name for name in archive.namelist() if name.endswith(".json")]
        target = next((name for name in candidates if name in names), None)
        if target is None:
            if len(names) == 1:
                target = names[0]
            else:
                sys.exit(f"Tree {tree_id} is not in {path}. Members: {', '.join(names[:10])}")
        return json.loads(archive.read(target))


def population_record(populations: dict, population_id: str) -> dict:
    if population_id in populations:
        return populations[population_id]
    if str(population_id) in populations:
        return populations[str(population_id)]
    sys.exit(f"Population {population_id} is missing from the selected PhyloWGS tree.")


def cellular_prevalence(population: dict) -> str:
    value = population.get("cellular_prevalence")
    if isinstance(value, list):
        if not value:
            return ""
        value = value[0]
    if value is None:
        return ""
    return f"{float(value):.8g}"


def configured_phylowgs_mode(root: Path) -> str:
    completed = subprocess.run(
        [sys.executable, str(root / "src" / "read_params.py"), "--format", "shell"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    for line in completed.stdout.splitlines():
        if line.startswith("PHYLOWGS_MODE="):
            return line.split("=", 1)[1].strip("'")
    return "smoke"


def phylowgs_dir(root: Path) -> tuple[Path, str]:
    base = root / "results" / "phylowgs"
    mode = configured_phylowgs_mode(root)
    preferred = base / mode
    if (preferred / "mutass.json").exists() and (preferred / "summ.json").exists():
        return preferred, mode
    for candidate in ("full", "smoke"):
        path = base / candidate
        if (path / "mutass.json").exists() and (path / "summ.json").exists():
            return path, candidate
    sys.exit("No PhyloWGS summ.json and mutass.json were found under results/phylowgs/.")


def extract_phylowgs(root: Path, mutations: list[dict[str, str]], sample: str) -> None:
    result_dir, mode = phylowgs_dir(root)
    summary = load_summary(result_dir / "summ.json")
    tree_id = choose_tree(summary)
    trees = summary["trees"]
    if tree_id not in trees:
        sys.exit(f"Selected PhyloWGS tree {tree_id} is not in summ.json.")
    mutass = load_mutass(result_dir / "mutass.json", tree_id)
    populations = trees[tree_id]["populations"]
    assignments = mutass["mut_assignments"]
    sid_to_pop: dict[str, str] = {}
    for population_id, content in assignments.items():
        for sid in content.get("ssms", []):
            sid_to_pop[sid] = str(population_id)

    ssm_map = {
        row["sid"]: row["mutation_id"]
        for row in read_tsv(root / "data" / "prepared" / "inputs" / "phylowgs" / "ssm_map.tsv")
    }
    pop_to_mutation: dict[str, tuple[str, str]] = {}
    for sid, population_id in sid_to_pop.items():
        mutation_id = ssm_map.get(sid)
        if mutation_id is None:
            sys.exit(f"PhyloWGS SSM {sid} is not in ssm_map.tsv")
        prevalence = cellular_prevalence(population_record(populations, population_id))
        pop_to_mutation[mutation_id] = (population_id, prevalence)

    rows = base_rows(mutations, "phylowgs")
    input_ids = set(ssm_map.values())
    for row in rows:
        found = pop_to_mutation.get(row["mutation_id"])
        if found is None:
            row["note"] = f"mcmc={mode};unassigned" if row["mutation_id"] in input_ids else "not in tool input"
            continue
        row["included"] = "true"
        row["cluster"] = found[0]
        row["ccf"] = found[1]
        row["note"] = f"mcmc={mode};tree={tree_id}"
    write_assignments(root / "results" / "assignments" / f"{sample}_phylowgs.tsv", rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool", choices=TOOLS)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    root = args.root.resolve()
    sample = sample_id(root)
    mutations = read_tsv(root / "data" / "prepared" / "mutations.tsv")
    selected = (args.tool,) if args.tool else TOOLS
    for tool in selected:
        if tool in ("decode", "mobster", "sciclone"):
            extract_r_tool(root, tool, mutations, sample)
        elif tool == "pyclonevi":
            extract_pyclone(root, mutations, sample)
        else:
            extract_phylowgs(root, mutations, sample)
        print(f"Wrote results/assignments/{sample}_{tool}.tsv")


if __name__ == "__main__":
    main()
