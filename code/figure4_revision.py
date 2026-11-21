#!/usr/bin/env python3
"""Final Figure 4: Python port of BFAST_publication.R Section 9.2 (62-gene OncoPrint at end of treatment).

Analysis copied from R, unchanged: anti_join(TRTDISC, SCRN, by = "pt_mutID"); variant_QC_pass; impact in {known, likely};
bep_62gene; Alterations = alteration_type_coding_type; distinct(PT, gene, Alterations); empty rows/columns removed.
No assay matching, CHIP filter, gene-level exclusion or manuscript-matching rule.

Reviewer 1 presentation fixes only: amplification yellow (#EDC948), rearrangement red (#FF0000) as in R; the legend is
fully visible; legend title spelled "Alterations". No loss category (Cohort A has no copy-number losses).

Usage: python figure4_revision.py [--data-dir DIR] [--outdir DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

DEFAULT_DATA_DIR = Path(__file__).resolve().parent / ".." / ".." / "data" / "BFAST-A"
DEFAULT_OAPL = "bfast_oapl_cohortA.withDBIannotations.clean.BlC3Td.csv"

# R `col` / `alter_fun`: key -> (colour, tile height fraction). Order = R `at`/legend order = drawing order.
TYPES = {
    "CN_amplification": ("#EDC948", 1.0),
    "RE_rearrangement": ("#FF0000", 1.0),
    "SV_frameshift": ("#008000", 0.33),
    "SV_missense": ("#0000FF", 0.33),
    "SV_nonframeshift": ("#79706E", 1.0),
    "SV_nonsense": ("#9D7660", 1.0),
    "SV_nonstop": ("#91DCEA", 1.0),
    "SV_promoter": ("#E39802", 1.0),
    "SV_splice": ("#BFBB60", 1.0),
}
BACKGROUND = "#CCCCCC"
LEGEND_TITLE = "Alterations"  # the R code has the typo "Alternations"


# ---------------------------------------------------------------------------
# Data: R Sections 9.1 / 9.2
# ---------------------------------------------------------------------------
def read_oapl(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False, na_values=["NA"])


def r_matrix(oapl: pd.DataFrame, only_62gene: bool) -> dict:
    """Return {'genes', 'patients', 'cells'}; cells[(gene, PT)] = ordered list of Alterations (pivot_wider list)."""
    scrn = oapl[oapl["VISIT"] == "SCRN"]
    trtdisc = oapl[oapl["VISIT"] == "TRTDISC"]
    tmp = trtdisc[~trtdisc["pt_mutID"].isin(set(scrn["pt_mutID"]))]          # anti_join(by = "pt_mutID")
    keep = tmp["variant_QC_pass"].eq("TRUE") & tmp["somatic_status_functional_impact"].isin(["known", "likely"])
    if only_62gene:
        keep &= tmp["bep_62gene"].eq("TRUE")
    sel = tmp[keep].copy()
    sel["Alterations"] = sel["alteration_type"] + "_" + sel["coding_type"]
    n_records = len(sel)                                                     # variant records before distinct()
    sel = sel[["PT", "gene", "Alterations"]].drop_duplicates()               # distinct()
    cells: dict[tuple[str, str], list[str]] = {}
    for pt, gene, alt in zip(sel["PT"], sel["gene"], sel["Alterations"]):
        cells.setdefault((gene, pt), []).append(alt)
    return {"genes": list(dict.fromkeys(sel["gene"])), "patients": list(dict.fromkeys(sel["PT"])), "cells": cells,
            "n_records": n_records}


def oncoprint_order(mat: dict) -> tuple[list[str], list[str]]:
    """ComplexHeatmap oncoPrint ordering: genes by decreasing number of alterations (stable); patients memo-sorted."""
    cells, rank = mat["cells"], {k: i for i, k in enumerate(TYPES)}
    count = {g: sum(len(cells.get((g, p), [])) for p in mat["patients"]) for g in mat["genes"]}
    genes = sorted(mat["genes"], key=lambda g: -count[g])

    def score(p):
        return [(1, -min(rank[a] for a in cells[(g, p)])) if (g, p) in cells else (0, 0) for g in genes]

    return genes, sorted(mat["patients"], key=score, reverse=True)


def matrix_fingerprint(mat: dict) -> str:
    rows = sorted(f"{p}|{g}|{';'.join(sorted(a))}" for (g, p), a in mat["cells"].items())
    return hashlib.sha256("\n".join(rows).encode()).hexdigest()[:16]


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def plot(mat: dict, genes: list[str], patients: list[str], outbase: Path, formats=("png", "pdf", "svg")) -> dict:
    n_pt, n_gene, cells = len(patients), len(genes), mat["cells"]
    cw, ch, gap = 0.36, 0.30, 0.04
    pct_w, names_w, bar_w, top_h, margin = 0.55, 0.95, 1.25, 0.95, 0.15
    legend_gap, key_w, row_h, text_size = 0.35, 0.22, 0.27, 8.5
    main_w, main_h = n_pt * cw, n_gene * ch
    keys = list(TYPES)

    # legend width (inches), measured from the longest label with the real renderer
    probe = plt.figure(figsize=(1, 1))
    renderer = probe.canvas.get_renderer()
    widest = max(probe.text(0, 0, k, fontsize=text_size).get_window_extent(renderer).width for k in keys) / probe.dpi
    plt.close(probe)
    legend_w = key_w + 0.08 + widest
    full_w = margin + pct_w + main_w + names_w + bar_w + legend_gap + legend_w + margin
    fig_w = full_w  # canvas includes the full legend width, so no label can be cut off
    bottom = 0.55  # room for the bar-chart axis labels
    fig_h = bottom + top_h + 0.1 + max(main_h, (len(keys) + 2) * row_h) + margin
    fig = plt.figure(figsize=(fig_w, fig_h))

    def axes_at(x, y, w, h):
        return fig.add_axes([x / fig_w, y / fig_h, w / fig_w, h / fig_h])

    body_x, body_y = margin + pct_w, fig_h - margin - top_h - 0.1 - main_h
    ax = axes_at(body_x, body_y, main_w, main_h)
    ax_top = axes_at(body_x, body_y + main_h + 0.1, main_w, top_h)
    ax_bar = axes_at(body_x + main_w + names_w, body_y, bar_w, main_h)

    for i, g in enumerate(genes):
        for j, p in enumerate(patients):
            ax.add_patch(Rectangle((j + gap, i + gap), 1 - 2 * gap, 1 - 2 * gap, fc=BACKGROUND, ec="none", lw=0))
            for k in keys:                                                  # alter_fun order: later overlays earlier
                if k in cells.get((g, p), []):
                    colour, frac = TYPES[k]
                    h = (1 - 2 * gap) * frac
                    ax.add_patch(Rectangle((j + gap, i + gap + ((1 - 2 * gap) - h) / 2), 1 - 2 * gap, h,
                                           fc=colour, ec="none", lw=0))
    ax.set_xlim(0, n_pt)
    ax.set_ylim(n_gene, 0)
    ax.set_xticks([])
    ax.set_yticks([i + 0.5 for i in range(n_gene)])
    ax.yaxis.tick_right()
    ax.set_yticklabels(genes, fontsize=9)
    ax.tick_params(length=0, pad=6)
    for s in ax.spines.values():
        s.set_visible(False)
    for i, g in enumerate(genes):                                           # % of the patients shown (ComplexHeatmap)
        pct = round(100 * sum(1 for p in patients if (g, p) in cells) / n_pt)
        ax.text(-0.15, i + 0.5, f"{pct}%", ha="right", va="center", fontsize=8.5)

    tmax = 1
    for j, p in enumerate(patients):                                         # alterations per patient
        b = 0
        for k in keys:
            c = sum(1 for g in genes if k in cells.get((g, p), []))
            if c:
                ax_top.bar(j + 0.5, c, bottom=b, width=1 - 2 * gap, color=TYPES[k][0], lw=0)
                b += c
        tmax = max(tmax, b)
    ax_top.set_xlim(0, n_pt)
    ax_top.set_ylim(0, tmax)
    ax_top.set_xticks([])
    ax_top.yaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=4))
    ax_top.tick_params(labelsize=8)
    for s in ("top", "right", "bottom"):
        ax_top.spines[s].set_visible(False)

    rmax = 1
    for i, g in enumerate(genes):                                            # alterations per gene
        x = 0
        for k in keys:
            c = sum(1 for p in patients if k in cells.get((g, p), []))
            if c:
                ax_bar.barh(i + 0.5, c, left=x, height=1 - 2 * gap, color=TYPES[k][0], lw=0)
                x += c
        rmax = max(rmax, x)
    ax_bar.set_ylim(n_gene, 0)
    ax_bar.set_xlim(0, rmax)
    ax_bar.set_yticks([])
    ax_bar.xaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=4))
    ax_bar.tick_params(labelsize=8)
    for s in ("top", "right", "left"):
        ax_bar.spines[s].set_visible(False)

    # legend as in R: title "Alternations", all nine entries, glyphs as drawn in the plot
    lx = margin + pct_w + main_w + names_w + bar_w + legend_gap
    ly = body_y + main_h
    fig.text(lx / fig_w, ly / fig_h, LEGEND_TITLE, fontsize=9, fontweight="bold", va="bottom")
    leg = axes_at(lx, ly - len(keys) * row_h - 0.05, legend_w, len(keys) * row_h)
    leg.set_xlim(0, legend_w)
    leg.set_ylim(0, len(keys) * row_h)
    leg.axis("off")
    for n, k in enumerate(keys):
        y = (len(keys) - 1 - n) * row_h
        leg.add_patch(Rectangle((0, y + 0.03), key_w, row_h - 0.06, fc=BACKGROUND, ec="none"))
        colour, frac = TYPES[k]
        hh = (row_h - 0.06) * frac
        leg.add_patch(Rectangle((0, y + 0.03 + (row_h - 0.06 - hh) / 2), key_w, hh, fc=colour, ec="none"))
        leg.text(key_w + 0.08, y + row_h / 2, k, fontsize=text_size, va="center", clip_on=False)

    fig.canvas.draw()
    clipped = clipped_text(fig)
    drawn = drawn_matrix(ax, genes, patients)
    for fmt in formats:
        fig.savefig(outbase.with_suffix("." + fmt), dpi=600 if fmt == "png" else None, facecolor="white")
    plt.close(fig)
    return {"clipped": clipped, "drawn": drawn, "size_in": (round(float(fig_w), 2), round(float(fig_h), 2))}


def clipped_text(fig) -> list[str]:
    r, c = fig.canvas.get_renderer(), fig.bbox
    out = []
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and t.get_text().strip():
            bb = t.get_window_extent(r)
            if bb.x0 < -0.5 or bb.y0 < -0.5 or bb.x1 > c.x1 + 0.5 or bb.y1 > c.y1 + 0.5:
                out.append(t.get_text())
    return out


def drawn_matrix(ax, genes, patients) -> set[tuple[str, str, str]]:
    """Reconstruct (gene, patient, alteration) from the rectangles actually drawn in the heatmap."""
    by_colour = {c.upper(): k for k, (c, _) in TYPES.items()}
    out = set()
    for p in ax.patches:
        colour = matplotlib.colors.to_hex(p.get_facecolor()).upper()
        if colour != BACKGROUND and colour in by_colour:
            out.add((genes[int(p.get_y())], patients[int(p.get_x())], by_colour[colour]))
    return out


# ---------------------------------------------------------------------------
EXPECTED = {"patients": 15, "genes": 9, "records": 21, "cells": 20}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    ap.add_argument("--oapl", default=DEFAULT_OAPL)
    ap.add_argument("--outdir", type=Path, default=Path("figure4_output"))
    args = ap.parse_args(argv)
    oapl = read_oapl(args.data_dir / args.oapl)
    args.outdir.mkdir(parents=True, exist_ok=True)

    mat = r_matrix(oapl, only_62gene=True)                                   # R Section 9.2
    genes, patients = oncoprint_order(mat)
    expected = {(g, p, a) for (g, p), alts in mat["cells"].items() for a in alts}
    res = plot(mat, genes, patients, args.outdir / "figure4")
    got = {"patients": len(patients), "genes": len(genes), "records": mat["n_records"], "cells": len(expected)}
    checks = {
        "counts equal the expected R result (15 patients, 9 genes, 21 variant records, 20 cells)": got == EXPECTED,
        "drawn matrix equals the R-derived matrix": res["drawn"] == expected,
        "legend and all text fully inside the canvas": not res["clipped"],
    }
    print(f"R Section 9.2: {got}; matrix fingerprint {matrix_fingerprint(mat)}; canvas {res['size_in']} in")
    print(f"gene order: {', '.join(genes)}")
    for name, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    if res["clipped"]:
        print("  clipped text:", res["clipped"])
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
