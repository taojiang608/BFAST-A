BFAST Cohort A -- Reviewer #3 Comment (3): Pearson correlation with 95% CI
==========================================================================

DELIVERABLES
- BFAST_ctDNA_SLD_correlation.py : Reproducible Python 3 analysis
- BFAST_ctDNA_SLD_requirements.txt : Package requirements
- BFAST_ctDNA_SLD_results/results.json : Machine-readable output
- BFAST_ctDNA_SLD_results/results_summary.txt : Human-readable output

PROVENANCE
The source R script BFAST_publication.R uses OSPL VISIT='SCRN' for baseline
ctDNA tumor fraction (cTF), log(cTF), and ADTR records flagged ABLFL='Y'
with PARAMCD='SUMLDIAM' for IRF SLD or PARAMCD='SD' for investigator SLD.
These are *different estimands*. The manuscript must identify which SLD
source was used for the reported r=0.24, p=0.024.

RUN
  python -m pip install -r BFAST_ctDNA_SLD_requirements.txt
  python BFAST_ctDNA_SLD_correlation.py

When ADTR is obtained from the authorized clinical ADaM data source:
  python BFAST_ctDNA_SLD_correlation.py --adtr adtr.csv --param both --make-plot

The script computes Pearson r, 2-sided p, exact paired n, Fisher-z 95% CI,
r^2, sensitivity Spearman rho, and plots for each baseline SLD definition.
It checks for repeated baseline clinical records rather than silently
averaging or duplicating patients, and exports aggregate results only.

STATISTICAL DEFINITION
  z = atanh(r)
  SE(z) = 1/sqrt(n-3)
  CI95(r) = tanh(z +/- 1.9599639845*SE(z))
  r^2 is variance explained for a simple unadjusted OLS fit with intercept.

VALIDATION NOTES
- The sample files have 87 positive baseline cTF measurements.
- The rounded r=0.24 and assumption n=87 yield CI approximately 0.031-0.429.
- p implied by r exactly 0.24 and n=87 is 0.02515, slightly different
  from the published p=0.024 because r has been rounded. That difference
  does not by itself indicate an error.
- The uploaded baseline cTF range is 0.06%-73.35%, whereas the manuscript
  currently reports 0.06%-32.90%. This discrepancy may reflect a different
  analysis subset or an exclusion; verify against the source ADTR merge
  before drawing conclusions about influence of possible outliers.

CONFIDENTIALITY
Patient-level rows, IDs and records are not written to analysis outputs.
Use only company-approved systems and study-specific data permissions for
any subsequent clinical database retrieval or analysis.
