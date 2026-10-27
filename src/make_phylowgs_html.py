#!/usr/bin/env python3
"""Build the one-page PhyloWGS HTML from a finished run.

Reads the summary, mutation list, and assignment archive written by
write_results.py, keeps at most 300 trees spread across posterior density,
and writes <sample>.html in the repository root. The page is SVG and CSS,
so it opens by double-click.

    python3 src/make_phylowgs_html.py
    python3 src/make_phylowgs_html.py results/phylowgs/full
    python3 src/make_phylowgs_html.py phylowgs_result/S8-noXY --out S8-noXY.html
"""

from __future__ import annotations

import argparse
import gzip
import html
import json
import re
import sys
import tempfile
import zipfile
from pathlib import Path

VAF_BIN = 0.05
PREV_BIN = 0.05
MIN_SSMS = 3
TOP_POPS = 3
EPSILON = 1e-6
POP_COLORS = ("#57534e", "#1f4e79", "#c47b2b", "#2f6f4e", "#8c3a4b", "#5b4b8a")


def load_json(path: Path):
    return json.loads(path.read_text())


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def histogram(values: list[float], width: float, start: float, stop: float) -> list[int]:
    count = max(1, int(round((stop - start) / width)))
    bins = [0] * count
    for value in values:
        if value < start or value > stop:
            continue
        index = min(count - 1, int((value - start) / width))
        bins[index] += 1
    return bins


def pop_color(pop_id: str) -> str:
    try:
        index = int(pop_id)
    except ValueError:
        index = 0
    return POP_COLORS[index % len(POP_COLORS)]


def vaf_of(ssm: dict) -> float | None:
    refs = ssm.get("ref_reads") or []
    totals = ssm.get("total_reads") or []
    if not refs or not totals or not totals[0]:
        return None
    return (totals[0] - refs[0]) / totals[0]


def cancer_population_count(populations: dict) -> int:
    count = 0
    for pop in populations.values():
        if pop.get("num_cnvs", 0) > 0 or pop.get("num_ssms", 0) >= MIN_SSMS:
            count += 1
    return count


def branch_ratio(tree: dict) -> float:
    branching = float(tree.get("branching_index") or 0)
    linearity = float(tree.get("linearity_index") or 0)
    return branching / (branching + linearity + EPSILON)


def layout_tree(structure: dict, populations: dict) -> list[dict]:
    children = {str(parent): [str(child) for child in kids] for parent, kids in structure.items()}
    ids = [str(pop_id) for pop_id in populations]
    child_ids = {child for kids in children.values() for child in kids}
    roots = [pop_id for pop_id in ids if pop_id not in child_ids] or ids[:1]
    rows: list[dict] = []
    cursor = {"x": 0}

    def place(node: str, depth: int) -> float:
        kids = children.get(node, [])
        if not kids:
            x = cursor["x"]
            cursor["x"] += 1
        else:
            x = mean([place(kid, depth + 1) for kid in kids])
        rows.append({"id": node, "x": x, "depth": depth})
        return x

    for root in roots:
        place(root, 0)
    return rows


def svg_histogram(
    bins: list[int],
    start: float,
    width: float,
    title: str,
    x_label: str,
    y_label: str,
    color: str,
) -> str:
    chart_w, chart_h = 820, 280
    left, right, top, bottom = 58, 16, 36, 46
    plot_w = chart_w - left - right
    plot_h = chart_h - top - bottom
    peak = max(bins) if bins else 1
    peak = max(peak, 1)
    bar_w = plot_w / max(len(bins), 1)
    parts = [
        f'<svg viewBox="0 0 {chart_w} {chart_h}" role="img">',
        f"<title>{html.escape(title)}</title>",
        f'<text x="{left}" y="22" class="chart-title">{html.escape(title)}</text>',
    ]
    for tick in range(5):
        value = peak * tick / 4
        y = top + plot_h - (plot_h * tick / 4)
        parts.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" class="grid"/>'
        )
        parts.append(
            f'<text x="{left - 8}" y="{y + 4:.1f}" class="tick" text-anchor="end">{int(round(value))}</text>'
        )
    for index, count in enumerate(bins):
        bar_h = plot_h * count / peak
        x = left + index * bar_w
        y = top + plot_h - bar_h
        parts.append(
            f'<rect x="{x + 1:.1f}" y="{y:.1f}" width="{max(bar_w - 2, 1):.1f}" '
            f'height="{bar_h:.1f}" fill="{color}"/>'
        )
    axis_y = top + plot_h
    parts.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{axis_y}" class="axis"/>')
    parts.append(
        f'<line x1="{left}" y1="{axis_y}" x2="{left + plot_w}" y2="{axis_y}" class="axis"/>'
    )
    stop = start + width * len(bins)
    for fraction in (0, 0.5, 1):
        label = start + (stop - start) * fraction
        x = left + plot_w * fraction
        parts.append(
            f'<text x="{x:.1f}" y="{axis_y + 18}" class="tick" text-anchor="middle">{label:.2f}</text>'
        )
    parts.append(
        f'<text x="{left + plot_w / 2:.1f}" y="{chart_h - 6}" class="axis-label" text-anchor="middle">'
        f"{html.escape(x_label)}</text>"
    )
    parts.append(
        f'<text x="16" y="{top + plot_h / 2:.1f}" class="axis-label" text-anchor="middle" '
        f'transform="rotate(-90 16 {top + plot_h / 2:.1f})">{html.escape(y_label)}</text>'
    )
    parts.append("</svg>")
    return "\n".join(parts)


def svg_columns(pairs: list[tuple[str, int]], title: str, x_label: str, y_label: str) -> str:
    chart_w, chart_h = 820, 280
    left, right, top, bottom = 58, 16, 36, 46
    plot_w = chart_w - left - right
    plot_h = chart_h - top - bottom
    peak = max((count for _, count in pairs), default=1)
    peak = max(peak, 1)
    slot = plot_w / max(len(pairs), 1)
    bar_w = min(slot * 0.64, 72)
    parts = [
        f'<svg viewBox="0 0 {chart_w} {chart_h}" role="img">',
        f"<title>{html.escape(title)}</title>",
        f'<text x="{left}" y="22" class="chart-title">{html.escape(title)}</text>',
    ]
    for tick in range(5):
        value = peak * tick / 4
        y = top + plot_h - (plot_h * tick / 4)
        parts.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" class="grid"/>'
        )
        parts.append(
            f'<text x="{left - 8}" y="{y + 4:.1f}" class="tick" text-anchor="end">{int(round(value))}</text>'
        )
    for index, (label, count) in enumerate(pairs):
        bar_h = plot_h * count / peak
        x = left + index * slot
        y = top + plot_h - bar_h
        parts.append(
            f'<rect x="{x + (slot - bar_w) / 2:.1f}" y="{y:.1f}" width="{bar_w:.1f}" '
            f'height="{bar_h:.1f}" fill="#1f4e79"/>'
        )
        parts.append(
            f'<text x="{x + slot / 2:.1f}" y="{top + plot_h + 18}" class="tick" text-anchor="middle">'
            f"{html.escape(label)}</text>"
        )
    axis_y = top + plot_h
    parts.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{axis_y}" class="axis"/>')
    parts.append(
        f'<line x1="{left}" y1="{axis_y}" x2="{left + plot_w}" y2="{axis_y}" class="axis"/>'
    )
    parts.append(
        f'<text x="{left + plot_w / 2:.1f}" y="{chart_h - 6}" class="axis-label" text-anchor="middle">'
        f"{html.escape(x_label)}</text>"
    )
    parts.append(
        f'<text x="16" y="{top + plot_h / 2:.1f}" class="axis-label" text-anchor="middle" '
        f'transform="rotate(-90 16 {top + plot_h / 2:.1f})">{html.escape(y_label)}</text>'
    )
    parts.append("</svg>")
    return "\n".join(parts)


def _padded(values: list[float]) -> tuple[float, float]:
    lo, hi = min(values), max(values)
    span = hi - lo or 0.05
    pad = span * 0.12
    return lo - pad, hi + pad


def svg_scatter(points: list[tuple[float, float, str]], best: str, title: str) -> str:
    chart_w, chart_h = 820, 460
    left, right, top, bottom = 64, 20, 36, 52
    plot_w = chart_w - left - right
    plot_h = chart_h - top - bottom
    xmin, xmax = _padded([point[0] for point in points])
    ymin, ymax = _padded([point[1] for point in points])

    def sx(value: float) -> float:
        return left + plot_w * (value - xmin) / (xmax - xmin)

    def sy(value: float) -> float:
        return top + plot_h * (1 - (value - ymin) / (ymax - ymin))

    def tick(value: float, span: float) -> str:
        digits = 2 if span >= 0.2 else 3
        return f"{value:.{digits}f}"

    parts = [
        f'<svg viewBox="0 0 {chart_w} {chart_h}" role="img">',
        f"<title>{html.escape(title)}</title>",
        f'<text x="{left}" y="22" class="chart-title">{html.escape(title)}</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" class="plot-bg"/>',
    ]
    for step in range(5):
        x_value = xmin + (xmax - xmin) * step / 4
        y_value = ymin + (ymax - ymin) * step / 4
        x = sx(x_value)
        y = sy(y_value)
        parts.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + plot_h}" class="grid"/>')
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" class="grid"/>')
        parts.append(
            f'<text x="{x:.1f}" y="{top + plot_h + 18}" class="tick" text-anchor="middle">'
            f"{tick(x_value, xmax - xmin)}</text>"
        )
        parts.append(
            f'<text x="{left - 8}" y="{y + 4:.1f}" class="tick" text-anchor="end">'
            f"{tick(y_value, ymax - ymin)}</text>"
        )
    for x_value, y_value, tree_id in points:
        if tree_id == best:
            continue
        parts.append(f'<circle cx="{sx(x_value):.1f}" cy="{sy(y_value):.1f}" r="3.2" class="point"/>')
    for x_value, y_value, tree_id in points:
        if tree_id != best:
            continue
        x, y = sx(x_value), sy(y_value)
        label_x = x + 12
        anchor = "start"
        if x > left + plot_w - 100:
            label_x = x - 12
            anchor = "end"
        parts.append(f'<line x1="{x - 9:.1f}" y1="{y:.1f}" x2="{x + 9:.1f}" y2="{y:.1f}" class="cross"/>')
        parts.append(f'<line x1="{x:.1f}" y1="{y - 9:.1f}" x2="{x:.1f}" y2="{y + 9:.1f}" class="cross"/>')
        parts.append(
            f'<text x="{label_x:.1f}" y="{y - 8:.1f}" class="best-label" text-anchor="{anchor}">'
            f"tree {html.escape(best)}</text>"
        )
    parts.append(
        f'<text x="{left + plot_w / 2:.1f}" y="{chart_h - 8}" class="axis-label" text-anchor="middle">'
        "Clustering index</text>"
    )
    parts.append(
        f'<text x="18" y="{top + plot_h / 2:.1f}" class="axis-label" text-anchor="middle" '
        f'transform="rotate(-90 18 {top + plot_h / 2:.1f})">BI / (BI + LI)</text>'
    )
    parts.append("</svg>")
    return "\n".join(parts)


def svg_phylogeny(tree: dict) -> str:
    structure = tree.get("structure") or {}
    populations = tree["populations"]
    nodes = layout_tree(structure, populations)
    by_id = {node["id"]: node for node in nodes}
    children = {str(parent): [str(child) for child in kids] for parent, kids in structure.items()}
    max_x = max(node["x"] for node in nodes)
    max_depth = max(node["depth"] for node in nodes)
    chart_w = max(820, int(220 * max_x + 280))
    chart_h = 120 + 130 * (max_depth + 1)
    top = 28
    offset = (chart_w - 220 * max_x) / 2

    def px(node: dict) -> float:
        return offset + node["x"] * 220

    def py(node: dict) -> float:
        return top + node["depth"] * 130 + 28

    parts = [
        f'<svg viewBox="0 0 {chart_w} {chart_h}" role="img">',
        "<title>Posterior mode tree</title>",
    ]
    for parent, kids in children.items():
        if parent not in by_id:
            continue
        for kid in kids:
            if kid not in by_id:
                continue
            parts.append(
                f'<line x1="{px(by_id[parent]):.1f}" y1="{py(by_id[parent]) + 28:.1f}" '
                f'x2="{px(by_id[kid]):.1f}" y2="{py(by_id[kid]) - 28:.1f}" class="edge"/>'
            )
    for node in nodes:
        pop = populations[node["id"]]
        prevalence = mean(pop.get("cellular_prevalence") or [0])
        x, y = px(node), py(node)
        title = "normal" if node["id"] == "0" else f"pop {node['id']}"
        detail = f"{prevalence:.3f}"
        counts = f"{pop.get('num_ssms', 0)} SSM, {pop.get('num_cnvs', 0)} CNV"
        parts.append(
            f'<rect x="{x - 62:.1f}" y="{y - 30:.1f}" width="124" height="60" rx="8" '
            f'fill="{pop_color(node["id"])}"/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{y - 8:.1f}" class="node-title" text-anchor="middle">{html.escape(title)}</text>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{y + 10:.1f}" class="node-detail" text-anchor="middle">{detail}</text>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{y + 24:.1f}" class="node-detail" text-anchor="middle">{html.escape(counts)}</text>'
        )
    parts.append("</svg>")
    return "\n".join(parts)


def coverage_note(view: Path) -> str:
    readme = view / "README.md"
    if not readme.exists():
        return ""
    text = readme.read_text()
    match = re.search(r"keeps (\d+) of (\d+) sampled trees", text)
    if match:
        return (
            f"The plots below use {match.group(1)} of {match.group(2)} sampled trees, "
            "spread across posterior density."
        )
    match = re.search(r"keeps all (\d+) sampled trees", text)
    if match:
        return f"The plots below use all {match.group(1)} sampled trees."
    return ""


def render(view: Path, meta_lead: str | None = None) -> str:
    sample_dirs = [path for path in (view / "data").iterdir() if path.is_dir()]
    if len(sample_dirs) != 1:
        sys.exit(f"Expected one sample directory in {view / 'data'}")
    sample_dir = sample_dirs[0]
    sample_id = sample_dir.name
    summary = load_json(sample_dir / f"{sample_id}.summ.json")
    muts = load_json(sample_dir / f"{sample_id}.muts.json")
    best = (view / "best_tree.txt").read_text().strip()
    assignment = load_json(sample_dir / f"{sample_id}.mutass" / f"{best}.json")
    tree = summary["trees"][best]
    ssms = muts.get("ssms") or {}

    all_vafs = [value for value in (vaf_of(ssm) for ssm in ssms.values()) if value is not None]
    vaf_bins = histogram(all_vafs, VAF_BIN, 0, 1)

    points = []
    pop_counts = []
    top_prevs = [[] for _ in range(TOP_POPS)]
    for tree_id, item in summary["trees"].items():
        points.append((float(item.get("clustering_index") or 0), branch_ratio(item), str(tree_id)))
        pop_counts.append(cancer_population_count(item["populations"]))
        ranked = sorted(
            item["populations"].items(),
            key=lambda pair: mean(pair[1].get("cellular_prevalence") or [0]),
            reverse=True,
        )
        cancer = ranked[1 : 1 + TOP_POPS]
        for index in range(TOP_POPS):
            if index < len(cancer):
                top_prevs[index].append(mean(cancer[index][1].get("cellular_prevalence") or [0]))

    count_min = min(pop_counts) if pop_counts else 0
    count_max = max(pop_counts) if pop_counts else 0
    count_pairs = []
    for value in range(count_min, count_max + 1):
        count_pairs.append((str(value), pop_counts.count(value)))

    assigned = assignment.get("mut_assignments") or {}
    pop_vaf_sections = []
    for pop_id in sorted(assigned, key=lambda item: int(item)):
        ids = assigned[pop_id].get("ssms") or []
        values = []
        for ssm_id in ids:
            ssm = ssms.get(ssm_id)
            if not ssm:
                continue
            value = vaf_of(ssm)
            if value is not None:
                values.append(value)
        if not values:
            pop_vaf_sections.append(
                f"<section><h2>Population {html.escape(pop_id)} VAF</h2>"
                "<p class='note'>This population has no SSMs.</p></section>"
            )
            continue
        bins = histogram(values, VAF_BIN, 0, 1)
        chart = svg_histogram(
            bins,
            0,
            VAF_BIN,
            f"Population {pop_id} VAF ({len(values)} SSMs)",
            "VAF",
            "SSMs",
            pop_color(pop_id),
        )
        pop_vaf_sections.append(f"<section><h2>Population {html.escape(pop_id)}</h2>{chart}</section>")

    prev_sections = []
    for index, values in enumerate(top_prevs, start=1):
        bins = histogram(values, PREV_BIN, 0, 1)
        chart = svg_histogram(
            bins,
            0,
            PREV_BIN,
            f"Cellular prevalence, cancer population {index} ({len(values)} trees)",
            "Cellular prevalence",
            "Trees",
            "#1f4e79",
        )
        prev_sections.append(chart)

    note = coverage_note(view)
    if meta_lead is None:
        meta_lead = (
            "Four chains, burn-in 1000, 2500 MCMC samples, "
            "5000 Metropolis–Hastings iterations, seeds 1001–1004."
        )
    prevalence = {
        pop_id: mean(pop.get("cellular_prevalence") or [0])
        for pop_id, pop in tree["populations"].items()
    }
    rows = []
    for pop_id in sorted(tree["populations"], key=lambda item: int(item)):
        pop = tree["populations"][pop_id]
        name = "normal" if pop_id == "0" else pop_id
        rows.append(
            "<tr>"
            f"<td>{html.escape(name)}</td>"
            f"<td>{prevalence[pop_id]:.3f}</td>"
            f"<td>{pop.get('num_ssms', 0)}</td>"
            f"<td>{pop.get('num_cnvs', 0)}</td>"
            "</tr>"
        )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>PhyloWGS {html.escape(sample_id)}</title>
<style>
:root {{ color-scheme: light; }}
body {{
  margin: 0;
  background: #f7f5f1;
  color: #1c1917;
  font: 16px/1.45 "Iowan Old Style", "Palatino Linotype", Palatino, Georgia, serif;
}}
main {{
  max-width: 920px;
  margin: 0 auto;
  padding: 32px 20px 64px;
}}
h1 {{ font-size: 2rem; line-height: 1.15; margin: 0 0 8px; }}
h2 {{ font-size: 1.25rem; margin: 0 0 8px; }}
p {{ margin: 0 0 12px; }}
.meta, .note {{ color: #44403c; }}
section {{
  background: #fff;
  border: 1px solid #e7e5e4;
  border-radius: 12px;
  padding: 18px 18px 8px;
  margin: 18px 0;
}}
svg {{ width: 100%; height: auto; display: block; }}
.chart-title, .node-title {{ font: 600 14px "Source Sans 3", "Segoe UI", sans-serif; fill: #1c1917; }}
.axis-label, .tick, .node-detail, .best-label {{
  font: 12px "Source Sans 3", "Segoe UI", sans-serif;
  fill: #44403c;
}}
.best-label {{ font-weight: 700; fill: #9a3412; }}
.grid {{ stroke: #e7e5e4; stroke-width: 1; }}
.axis {{ stroke: #1c1917; stroke-width: 1.25; }}
.plot-bg {{ fill: #fafaf9; }}
.point {{ fill: #1f4e79; opacity: 0.55; }}
.cross {{ stroke: #c2410c; stroke-width: 2.5; }}
.edge {{ stroke: #44403c; stroke-width: 1.5; }}
.node-title, .node-detail {{ fill: #fff; }}
table {{
  width: 100%;
  border-collapse: collapse;
  font: 15px "Source Sans 3", "Segoe UI", sans-serif;
  margin: 8px 0 16px;
}}
th, td {{ text-align: left; padding: 6px 8px; border-bottom: 1px solid #e7e5e4; }}
th {{ color: #44403c; font-weight: 650; }}
</style>
</head>
<body>
<main>
<h1>PhyloWGS {html.escape(sample_id)}</h1>
<p class="meta">{html.escape(meta_lead)} {html.escape(note)} The marked tree is the highest-density tree, {html.escape(best)}.</p>

<section>
<h2>Variant allele frequency</h2>
<p class="note">{len(all_vafs)} mutations.</p>
{svg_histogram(vaf_bins, 0, VAF_BIN, "All mutations", "VAF", "Mutations", "#1f4e79")}
</section>

<section>
<h2>Sampled trees</h2>
<p class="note">Each point is one retained tree. The cross is tree {html.escape(best)}.</p>
{svg_scatter(points, best, "Clustering index against branching")}
{svg_columns(count_pairs, f"Cancer populations ({len(pop_counts)} trees)", "Populations", "Trees")}
</section>

<section>
<h2>Cellular prevalence</h2>
<p class="note">For each tree, populations are ordered by cellular prevalence. The normal root is left out. These are the next three.</p>
{''.join(prev_sections)}
</section>

<section>
<h2>Tree {html.escape(best)}</h2>
<p class="note">Population 0 is the normal root and has no SSMs. Labels are cellular prevalence, then SSM and CNV counts. The VAF plots below are the cancer populations.</p>
{svg_phylogeny(tree)}
<table>
<thead><tr><th>Population</th><th>Cellular prevalence</th><th>SSMs</th><th>CNVs</th></tr></thead>
<tbody>
{''.join(rows)}
</tbody>
</table>
</section>

{''.join(pop_vaf_sections)}
</main>
</body>
</html>
"""


ROOT = Path(__file__).resolve().parents[1]
MAX_TREES = 300
MAX_SUMMARY_BYTES = 30 * 1024 * 1024


def is_view(path: Path) -> bool:
    return (path / "best_tree.txt").is_file() and (path / "data").is_dir()


def find_summary(directory: Path) -> Path | None:
    if not directory.is_dir():
        return None
    found: list[Path] = []
    for pattern in ("*.summ.json.gz", "*.summ.json", "summ.json.gz", "summ.json"):
        found.extend(sorted(directory.glob(pattern)))
    for path in found:
        if path.is_file() and path.stat().st_size > 0:
            return path
    return None


def load_json_bytes(path: Path) -> tuple[dict, int]:
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return json.loads(raw), len(raw)


def companion(summary: Path, kind: str) -> Path:
    stem = summary.name
    prefix = ""
    for suffix in (".summ.json.gz", ".summ.json"):
        if stem.endswith(suffix):
            prefix = stem[: -len(suffix)]
            break
    if kind == "muts":
        names = [f"{prefix}.muts.json.gz", f"{prefix}.muts.json", "muts.json.gz", "muts.json", "mutlist.json"]
    else:
        names = [f"{prefix}.mutass.zip", f"{prefix}.mutass.json", "mutass.zip", "mutass.json"]
    seen: set[str] = set()
    for name in names:
        if not name or name in seen or name.startswith("."):
            continue
        seen.add(name)
        path = summary.parent / name
        if path.is_file() and path.stat().st_size > 0:
            return path
    sys.exit(f"Could not find the {kind} file next to {summary}")


def sample_id_of(summary: Path, params: str) -> str:
    for line in params.splitlines():
        if line.startswith("sample="):
            return line.split("=", 1)[1].strip()
    stem = summary.name
    for suffix in (".summ.json.gz", ".summ.json"):
        if stem.endswith(suffix):
            prefix = stem[: -len(suffix)]
            if prefix:
                return prefix
    return "S8-noXY"


def format_seeds(text: str) -> str:
    parts = [part.strip() for part in text.split(",") if part.strip()]
    try:
        numbers = [int(part) for part in parts]
    except ValueError:
        return ", ".join(parts) or text
    if len(numbers) >= 2 and numbers == list(range(numbers[0], numbers[0] + len(numbers))):
        return f"{numbers[0]}–{numbers[-1]}"
    return ", ".join(parts)


def meta_from_params(params: str) -> str | None:
    values = {}
    for line in params.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    if not {"chains", "burnin", "mcmc", "mh", "seeds"} <= values.keys():
        return None
    chain_count = values["chains"]
    chain_word = {"1": "One chain", "2": "Two chains", "3": "Three chains", "4": "Four chains"}.get(
        chain_count, f"{chain_count} chains"
    )
    return (
        f"{chain_word}, burn-in {values['burnin']}, {values['mcmc']} MCMC samples, "
        f"{values['mh']} Metropolis–Hastings iterations, seeds {format_seeds(values['seeds'])}."
    )


def thin_summary(summary: dict, full_bytes: int) -> tuple[dict, str, int, int, bool]:
    trees = summary.get("trees") or {}
    densities = summary.get("tree_densities") or {}
    if not trees:
        sys.exit("Summary JSON contains no trees.")

    def density(tree_id: str) -> float:
        if tree_id in densities:
            return float(densities[tree_id])
        return float("-inf")

    ranked = sorted(trees, key=density)
    best = ranked[-1]
    full_count = len(trees)
    kept = list(ranked)
    thinned = False
    if full_count > MAX_TREES or full_bytes > MAX_SUMMARY_BYTES:
        thinned = True
        if full_count <= MAX_TREES:
            chosen = ranked
        else:
            chosen = []
            last = full_count - 1
            for index in range(MAX_TREES):
                pick = round(index * last / (MAX_TREES - 1))
                chosen.append(ranked[pick])
            if best not in chosen:
                chosen.append(best)
        seen: set[str] = set()
        kept = []
        for tree_id in chosen:
            if tree_id not in seen:
                seen.add(tree_id)
                kept.append(tree_id)
        summary["trees"] = {tree_id: trees[tree_id] for tree_id in kept}
        summary["tree_densities"] = {
            tree_id: densities[tree_id] for tree_id in kept if tree_id in densities
        }
    return summary, best, len(kept), full_count, thinned


def stage_result(result_dir: Path) -> tuple[Path, str, tempfile.TemporaryDirectory[str]]:
    summary_path = find_summary(result_dir)
    if summary_path is None:
        sys.exit(f"No PhyloWGS summary JSON in {result_dir}")
    params_path = result_dir / "params.txt"
    params = params_path.read_text() if params_path.is_file() else ""
    sample_id = sample_id_of(summary_path, params)
    summary, full_bytes = load_json_bytes(summary_path)
    summary, best, kept, full_count, thinned = thin_summary(summary, full_bytes)
    muts_path = companion(summary_path, "muts")
    mutass_path = companion(summary_path, "mutass")

    temp = tempfile.TemporaryDirectory(prefix="phylowgs-html-")
    stage = Path(temp.name)
    data_dir = stage / "data" / sample_id
    mutass_dir = data_dir / f"{sample_id}.mutass"
    mutass_dir.mkdir(parents=True)
    with (data_dir / f"{sample_id}.summ.json").open("w") as handle:
        json.dump(summary, handle, separators=(",", ":"))
    muts_raw = muts_path.read_bytes()
    if muts_raw[:2] == b"\x1f\x8b":
        muts_raw = gzip.decompress(muts_raw)
    (data_dir / f"{sample_id}.muts.json").write_bytes(muts_raw)

    member = None
    with zipfile.ZipFile(mutass_path) as archive:
        target = f"{best}.json"
        for name in archive.namelist():
            base = name.replace("\\", "/").rstrip("/").split("/")[-1]
            if base == target:
                member = name
                break
        if member is None:
            sys.exit(f"Assignment for tree {best} is not in {mutass_path}")
        (mutass_dir / target).write_bytes(archive.read(member))
    (stage / "best_tree.txt").write_text(best + "\n")
    if thinned:
        sentence = (
            f"The summary keeps {kept} of {full_count} sampled trees, "
            "spread across posterior density, including the highest-density tree."
        )
    else:
        sentence = f"The summary keeps all {full_count} sampled trees."
    (stage / "README.md").write_text(sentence + "\n")
    print(f"best_tree={best} kept={kept} of={full_count} thinned={thinned}")
    return stage, meta_from_params(params) or "", temp


def discover(root: Path) -> Path:
    for path in (
        root / "results" / "phylowgs" / "full",
        root / "phylowgs_result" / "S8-noXY",
        root / "phylowgs_view",
    ):
        if is_view(path) or find_summary(path):
            return path
    sys.exit(
        "No PhyloWGS result found. Pass the directory that contains the summary JSON."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Write the one-page PhyloWGS HTML.")
    parser.add_argument(
        "directory",
        nargs="?",
        type=Path,
        help="Result directory from write_results.py, or an assembled phylowgs_view directory.",
    )
    parser.add_argument("--out", type=Path, help="Output HTML path. Defaults to <repo>/<sample>.html.")
    args = parser.parse_args()
    source = args.directory if args.directory is not None else discover(ROOT)
    source = source.resolve()
    temp = None
    meta_lead = None
    if is_view(source):
        view = source
    else:
        view, meta_lead, temp = stage_result(source)
    try:
        sample_dirs = [path.name for path in (view / "data").iterdir() if path.is_dir()]
        if len(sample_dirs) != 1:
            sys.exit(f"Expected one sample under {view / 'data'}")
        target = args.out if args.out is not None else ROOT / f"{sample_dirs[0]}.html"
        target.write_text(render(view, meta_lead or None))
    finally:
        if temp is not None:
            temp.cleanup()
    print(f"Wrote {target}")


if __name__ == "__main__":
    main()
