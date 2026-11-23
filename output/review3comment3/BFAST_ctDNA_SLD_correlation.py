#!/usr/bin/env python3
"""BFAST Cohort A: baseline log(cTF) versus baseline SLD.

Purpose
-------
Reproduce Reviewer #3, Comment (3): Pearson correlation coefficient, two-sided
p-value, 95% confidence interval using Fisher's z transformation, and r squared.

IMPORTANT DATA REQUIREMENT
--------------------------
The uploaded OSPL plasma sample CSV has baseline cTF for 87 unique patients.
The supplied ADaM ADRS and ADTTE files do not contain baseline SLD. The
original BFAST R code retrieves SLD from the clinical ADTR ADaM dataset:
  ABLFL == 'Y'; PARAMCD == 'SUMLDIAM' (IRF) OR 'SD' (investigator);
  baseline SLD = numeric AVALC; PT = final SUBJID component.

Until ADTR is supplied, this program produces ONLY a *provisional* Fisher-z
confidence interval using the manuscript's rounded r = 0.24 and assuming all
87 baseline patients have paired baseline SLD. It does NOT recreate r or p
from patient-level paired SLD measurements.

Privacy
-------
This program writes aggregate outputs only; patient-level datasets are never
exported. Run in an authorized secure environment under your applicable
clinical-study and company data-handling rules.

Examples
--------
python BFAST_ctDNA_SLD_correlation.py
python BFAST_ctDNA_SLD_correlation.py --adtr adtr.csv --param both
python BFAST_ctDNA_SLD_correlation.py --adtr adtr.sas7bdat --param SUMLDIAM
python BFAST_ctDNA_SLD_correlation.py --adtr adtr.csv --param SD --make-plot
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


REFERENCE_R = 0.24
REFERENCE_P = 0.024
CONFIDENCE_LEVEL = 0.95


def fisher_ci(r: float, n: int, confidence_level: float = 0.95) -> tuple[float, float]:
    """Two-sided Fisher-z CI for population Pearson correlation, matching R cor.test()."""
    if n <= 3 or abs(r) >= 1.0:
        raise ValueError("Fisher-z CI needs n > 3 and -1 < r < 1")
    alpha = 1.0 - confidence_level
    margin = stats.norm.ppf(1 - alpha / 2) / math.sqrt(n - 3)
    z = np.arctanh(r)
    return (float(np.tanh(z - margin)), float(np.tanh(z + margin)))


def read_dataset(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Dataset not found: {path}")
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, low_memory=False)
    if suffix in {".xls", ".xlsx"}:
        return pd.read_excel(path)
    if suffix == ".sas7bdat":
        # pandas.read_sas supports standard SAS7BDAT exports.
        return pd.read_sas(path, format="sas7bdat", encoding="utf-8")
    raise ValueError(f"Unsupported dataset format: {path.name}")


def normalize_pt(values: pd.Series) -> pd.Series:
    """Normalize patient identifier to the numerical suffix, consistent with the R script."""
    strings = values.astype("string").str.strip()
    suffix = strings.str.extract(r"(\d+)(?:\.0+)?$", expand=False)
    return pd.to_numeric(suffix, errors="coerce").astype("Int64")


def baseline_ctf(ospl: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    required = {"PT", "VISIT", "cTF"}
    missing = required - set(ospl.columns)
    if missing:
        raise ValueError(f"OSPL is missing required columns: {sorted(missing)}")

    scrn = ospl.loc[ospl["VISIT"].astype(str).str.upper().str.strip().eq("SCRN")].copy()
    scrn["PT"] = normalize_pt(scrn["PT"])
    scrn["cTF"] = pd.to_numeric(scrn["cTF"], errors="coerce")

    if scrn["PT"].isna().any():
        raise ValueError("Could not parse some baseline patient IDs from OSPL.")
    if scrn["PT"].duplicated().any():
        raise ValueError("Multiple SCRN rows found for a patient. Resolve assay/visit selection first.")

    audit = {
        "baseline_plasma_rows": int(len(scrn)),
        "unique_patients_with_screening_plasma": int(scrn["PT"].nunique()),
        "nonmissing_ctf": int(scrn["cTF"].notna().sum()),
        "positive_ctf": int(scrn["cTF"].gt(0).sum()),
        "zero_or_negative_ctf": int(scrn["cTF"].le(0).sum()),
        "missing_ctf": int(scrn["cTF"].isna().sum()),
    }
    valid = scrn.loc[scrn["cTF"].gt(0), ["PT", "cTF"]].copy()
    if valid.empty:
        raise ValueError("There are no positive baseline cTF values for a log transform.")
    valid["logcTF"] = np.log(valid["cTF"].to_numpy(dtype=float))
    return valid, audit


def baseline_sld(adtr: pd.DataFrame, paramcd: str) -> tuple[pd.DataFrame, dict]:
    required = {"SUBJID", "ABLFL", "PARAMCD"}
    missing = required - set(adtr.columns)
    if missing:
        raise ValueError(f"ADTR is missing required columns: {sorted(missing)}")
    if "AVALC" not in adtr.columns:
        raise ValueError("ADTR needs AVALC, as in the source BFAST R script.")

    selected = adtr.loc[
        adtr["ABLFL"].astype(str).str.strip().str.upper().eq("Y") &
        adtr["PARAMCD"].astype(str).str.strip().str.upper().eq(paramcd)
    ].copy()
    selected["PT"] = normalize_pt(selected["SUBJID"])
    selected["SLD"] = pd.to_numeric(selected["AVALC"], errors="coerce")
    if selected["PT"].isna().any():
        raise ValueError(f"Cannot parse all baseline {paramcd} ADTR SUBJIDs")
    if selected["PT"].duplicated().any():
        raise ValueError(
            f"ADTR has more than one baseline row per subject for {paramcd}. "
            "Check evaluator/assessment metadata and match source R definitions before aggregating."
        )
    audit = {
        "paramcd": paramcd,
        "baseline_ADTR_rows": int(len(selected)),
        "numeric_baseline_sld": int(selected["SLD"].notna().sum()),
    }
    return selected[["PT", "SLD"]], audit


def analyze_pairing(ctf: pd.DataFrame, adtr: pd.DataFrame, paramcd: str) -> tuple[dict, pd.DataFrame]:
    sld, sld_audit = baseline_sld(adtr, paramcd)
    matched = ctf.merge(sld, how="inner", on="PT", validate="one_to_one")
    paired = matched.dropna(subset=["SLD", "logcTF"]).copy()
    if len(paired) <= 3:
        raise ValueError(f"Only {len(paired)} complete baseline pairs for {paramcd}")
    if paired["SLD"].nunique() < 2 or paired["logcTF"].nunique() < 2:
        raise ValueError(f"At least one paired variable has zero variance for {paramcd}")

    x = paired["SLD"].to_numpy(dtype=float)
    y = paired["logcTF"].to_numpy(dtype=float)
    result = stats.pearsonr(x, y)
    r = float(result.statistic)
    n = len(paired)
    lo, hi = fisher_ci(r, n)
    reg = stats.linregress(x, y)
    sp = stats.spearmanr(x, y)
    results = {
        "status": "reproduced_from_paired_baseline_SLD_and_ctf",
        **sld_audit,
        "matched_subjects_before_numeric_SLD_filter": int(len(matched)),
        "complete_pairs_n": int(n),
        "pearson_r": r,
        "pearson_p_two_sided": float(result.pvalue),
        "pearson_r_95CI_fisher": [lo, hi],
        "r_squared": r * r,
        "variance_explained_percent": 100 * r * r,
        "ols_intercept": float(reg.intercept),
        "ols_slope_logctf_per_unit_sld": float(reg.slope),
        "spearman_rho_sensitivity": float(sp.statistic),
        "spearman_p_two_sided_sensitivity": float(sp.pvalue),
        "r_difference_from_manuscript_rounded_0_24": r - REFERENCE_R,
    }
    return results, paired


def make_plot(paired: pd.DataFrame, result: dict, output_path: Path) -> None:
    import matplotlib.pyplot as plt

    x = paired["SLD"].to_numpy(dtype=float)
    y = paired["logcTF"].to_numpy(dtype=float)
    slope = result["ols_slope_logctf_per_unit_sld"]
    intercept = result["ols_intercept"]
    xx = np.linspace(float(np.min(x)), float(np.max(x)), 100)

    fig, ax = plt.subplots(figsize=(7.0, 5.0), dpi=180)
    ax.scatter(x, y, s=32, alpha=0.75)
    ax.plot(xx, slope * xx + intercept, linewidth=1.75)
    ci = result["pearson_r_95CI_fisher"]
    ax.set_title(f"Baseline {result['paramcd']}: log(cTF) vs SLD")
    ax.set_xlabel("Baseline sum of longest diameters (SLD)")
    ax.set_ylabel("Natural log of baseline cTF")
    ax.text(
        0.04, 0.96,
        f"n={result['complete_pairs_n']}  r={result['pearson_r']:.3f} "
        f"(95% CI {ci[0]:.3f}, {ci[1]:.3f})\n"
        f"r²={result['r_squared']:.3f}  p={result['pearson_p_two_sided']:.3g}",
        transform=ax.transAxes, va="top", fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "0.8", "alpha": 0.9},
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def provisional_report(n_assumed: int) -> dict:
    lo, hi = fisher_ci(REFERENCE_R, n_assumed)
    t = REFERENCE_R * math.sqrt((n_assumed - 2) / (1 - REFERENCE_R ** 2))
    p_implied = 2 * stats.t.sf(abs(t), df=n_assumed - 2)
    return {
        "status": "PROVISIONAL_ONLY_NOT_REPRODUCED",
        "reason": "ADTR with baseline SLD is not supplied; only the rounded manuscript r is known.",
        "assumed_n_complete_pairs": int(n_assumed),
        "reported_r_rounded": REFERENCE_R,
        "reported_p_rounded": REFERENCE_P,
        "95CI_fisher_from_reported_r_and_assumed_n": [lo, hi],
        "r_squared_from_reported_rounded_r": REFERENCE_R ** 2,
        "variance_explained_percent_from_reported_r": 100 * REFERENCE_R ** 2,
        "p_implied_by_rounded_r_and_assumed_n": float(p_implied),
        "interpretation": "An unadjusted weak positive linear correlation; no causal inference.",
        "warning": "Confirm n, r, and CI against the original paired ADTR data before journal submission.",
    }


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ospl", type=Path, default=script_dir / "bfast_oSpl_cohortA.withDBIannotations.BlC3Td.csv")
    p.add_argument("--adtr", type=Path, default=None, help="ADTR ADaM clinical dataset with baseline SLD (CSV/XLSX/SAS7BDAT).")
    p.add_argument("--param", choices=["SUMLDIAM", "SD", "both"], default="both",
                   help="SUMLDIAM=IRF; SD=investigator. Source R code analyzes both.")
    p.add_argument("--output-dir", type=Path, default=script_dir / "BFAST_ctDNA_SLD_results")
    p.add_argument("--make-plot", action="store_true", help="Save scatterplot if ADTR is available.")
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    ospl = read_dataset(args.ospl)
    ctf, audit = baseline_ctf(ospl)
    findings: dict = {
        "study": "BFAST Cohort A - Reviewer 3 Comment 3",
        "methods": "Screening (SCRN) positive cTF, natural logarithm; ADTR ABLFL=Y; Pearson r; 95% CI by Fisher's z; unadjusted r-squared.",
        "input_ctf_file": args.ospl.name,
        "plasma_data_audit": audit,
        "analysis": {},
    }

    if args.adtr is None:
        findings["analysis"]["provisional"] = provisional_report(audit["positive_ctf"])
    else:
        adtr = read_dataset(args.adtr)
        findings["input_adtr_file"] = args.adtr.name
        for paramcd in (["SUMLDIAM", "SD"] if args.param == "both" else [args.param]):
            results, paired = analyze_pairing(ctf, adtr, paramcd)
            findings["analysis"][paramcd] = results
            if args.make_plot:
                make_plot(paired, results, args.output_dir / f"baseline_ctf_sld_{paramcd}.png")

    result_path = args.output_dir / "results.json"
    result_path.write_text(json.dumps(findings, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report_lines = [
        "BFAST Cohort A - Pearson correlation of screening log(cTF) and baseline SLD",
        "=" * 72,
        f"Screening plasma rows: {audit['baseline_plasma_rows']}",
        f"Unique screening patients: {audit['unique_patients_with_screening_plasma']}",
        f"Positive cTF values (eligible for ln): {audit['positive_ctf']}",
        "",
    ]
    for key, analysis in findings["analysis"].items():
        report_lines.append(f"[{key}] {analysis['status']}")
        if key == "provisional":
            low, high = analysis["95CI_fisher_from_reported_r_and_assumed_n"]
            report_lines += [
                f"Assumed paired n: {analysis['assumed_n_complete_pairs']}",
                f"Reported rounded r: {analysis['reported_r_rounded']:.2f}",
                f"Provisional 95% CI: ({low:.3f}, {high:.3f})",
                f"r^2 from rounded r: {analysis['r_squared_from_reported_rounded_r']:.4f}",
                f"Reported p: {analysis['reported_p_rounded']:.3f}",
                f"p implied by rounded r/n: {analysis['p_implied_by_rounded_r_and_assumed_n']:.4f}",
                "WARNING: NOT A REPRODUCTION WITHOUT ADTR / COMPLETE PAIRS",
            ]
        else:
            low, high = analysis["pearson_r_95CI_fisher"]
            report_lines += [
                f"Paired n: {analysis['complete_pairs_n']}",
                f"Pearson r: {analysis['pearson_r']:.6f}",
                f"95% CI: ({low:.6f}, {high:.6f})",
                f"p (two-sided): {analysis['pearson_p_two_sided']:.6g}",
                f"r^2: {analysis['r_squared']:.6f}",
                f"Spearman sensitivity rho: {analysis['spearman_rho_sensitivity']:.6f}",
            ]
        report_lines.append("")
    summary_path = args.output_dir / "results_summary.txt"
    summary_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print("\n".join(report_lines))
    print(f"Results saved to: {result_path} and {summary_path}")


if __name__ == "__main__":
    main()
