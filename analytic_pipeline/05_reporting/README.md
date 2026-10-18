# Manuscript reporting

## Available

- `30_build_tables.R` runs the archived reviewer-facing R implementation at `../implementation/30_build_SPINE_TRACT_final_tables_PUBLIC_REVIEWER_v4_1.R`.
- The table builder generates three landscape main tables, eTables 1–10 (including fold-specific analytic sample sizes), and machine-readable source CSVs.
- It requires the local **repaired diagnosis** cohort and frozen reporting artifacts. It does not run model training or retrieve ACS PUF records.

## Pending archive

- Final Figure 1–3 and eFigure 1–6 R source, after the final visual QA changes are incorporated.
- TRIPOD+AI subgroup bootstrap reporting source after final execution/QA.
- Full package/environment manifests and deidentified aggregate QC summaries.

Until the corresponding final source files are archived, figure or subgroup results in the manuscript should be traced to the authors' locked local analysis package, not claimed to be reproducible directly from this public repository.
