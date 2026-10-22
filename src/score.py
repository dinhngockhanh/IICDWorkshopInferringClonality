#!/usr/bin/env python3
"""Placeholder for the workshop score.

The scoring rules are not specified yet. This script only checks that the
assignment tables exist and then stops. It does not compute SMC-Het or any
other score.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    assignment_dir = ROOT / "results" / "assignments"
    tables = sorted(assignment_dir.glob("*.tsv")) if assignment_dir.exists() else []
    print("评分规则尚未给出，当前不计算分数。")
    if not tables:
        print("还没有 results/assignments/*.tsv。先运行各算法和 src/extract_results.py。")
        sys.exit(0)
    print("已找到这些统一结果表，可供以后评分：")
    for path in tables:
        print(f"  {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
