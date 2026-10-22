#!/usr/bin/env python3
"""Print config/params.toml as shell assignments."""

from __future__ import annotations

import argparse
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def shell_value(value) -> str:
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, str):
        return "'" + value.replace("'", "'\"'\"'") + "'"
    return str(value)


def assignments(config: dict) -> dict[str, object]:
    mode = config["phylowgs"]["mode"]
    if mode not in config["phylowgs"]:
        sys.exit(f"phylowgs.mode must be smoke or full, got {mode}")
    chain = config["phylowgs"][mode]
    return {
        "SAMPLE_ID": config["sample_id"],
        "SEED": config["seed"],
        "DECODE_NEUTRAL_TAIL": config["decode"]["neutral_tail"],
        "DECODE_COMPUTE_PARALLEL": config["decode"]["compute_parallel"],
        "MOBSTER_PARALLEL": config["mobster"]["parallel"],
        "MOBSTER_K_MIN": config["mobster"]["k_min"],
        "MOBSTER_K_MAX": config["mobster"]["k_max"],
        "SCICLONE_MINIMUM_DEPTH": config["sciclone"]["minimum_depth"],
        "SCICLONE_ALONG_MARGINS": config["sciclone"]["do_clustering_along_margins"],
        "PYCLONE_CLUSTERS": config["pyclonevi"]["clusters"],
        "PYCLONE_DENSITY": config["pyclonevi"]["density"],
        "PYCLONE_RESTARTS": config["pyclonevi"]["restarts"],
        "PHYLOWGS_MODE": mode,
        "PHYLOWGS_CHAINS": config["phylowgs"]["chains"],
        "PHYLOWGS_BURNIN": chain["burnin"],
        "PHYLOWGS_MCMC": chain["mcmc"],
        "PHYLOWGS_MH": chain["mh"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=["shell"], default="shell")
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "params.toml")
    args = parser.parse_args()
    with args.config.open("rb") as handle:
        config = tomllib.load(handle)
    for key, value in assignments(config).items():
        print(f"{key}={shell_value(value)}")


if __name__ == "__main__":
    main()
