# SEPTRALU re-treatment: reproducible analysis

De-identified data and R code to reproduce the main tables and figures of the
study of [177Lu]Lu-DOTATATE re-treatment (R-PRRT) in metastatic neuroendocrine
neoplasms within the SEPTRALU registry.

## Data

`septralu_retreatment.rds` contains one row per patient (n = 134). It is fully
de-identified: it holds no patient identifiers, no centre, and no dates — only
pre-computed intervals (in months) and analysis variables.

| Variable | Description |
|---|---|
| `age` | Age at first R-PRRT cycle (years) |
| `sex` | Sex |
| `primary_site` | Primary tumour site (Pancreas / Ileum / Lung / Other; `NA` = unknown primary) |
| `primary_site_pancreas` | Pancreas vs others (unknown grouped with others) |
| `ki67` | Ki-67 index (%) |
| `grade` | WHO grade |
| `ecog` | ECOG performance status (0–3) |
| `functioning_tumor` | Functioning tumour (No / Yes) |
| `peritoneal_mets`, `liver_mets` | Metastatic sites (No / Yes) |
| `pet_ga_heterogeneity` | Intratumoural 68Ga-PET uptake heterogeneity (No / Yes) |
| `retreatment_interval_months` | Interval from last I-PRRT cycle to first R-PRRT cycle |
| `metastasis_to_retreatment_months` | Interval from metastasis diagnosis to R-PRRT |
| `n_cycles` | Number of R-PRRT cycles (1–4) |
| `cycle2_month`, `cycle3_month`, `cycle4_month` | Time from the first R-PRRT cycle to each subsequent cycle (months) |
| `dose_reduced_c1` | Cycle 1 dose (Standard 7.4 GBq / Reduced) |
| `discontinuation_reason` | Reason for discontinuation |
| `recist` | RECIST 1.1 response (CR / PR / SD / PD / Non-evaluable) |
| `clinical_response` | Clinical response |
| `biochemical_response_cga`, `biochemical_response_secreted` | Biochemical response |
| `tox_nausea`, `tox_vomiting`, `tox_urticaria`, `tox_hematologic`, `tox_nephrotoxicity`, `tox_other` | Maximum CTCAE grade per patient (0–4) |
| `tox_overall_max` | Maximum toxicity grade across categories |
| `tox_evaluable` | Toxicity data available |
| `intercurrent_systemic_therapy` | Systemic antitumour therapy (chemotherapy/targeted/immunotherapy) started between the last I-PRRT and the first R-PRRT cycle (No / Yes) |
| `os_time`, `os_event` | Overall survival (months; 1 = death) |
| `pfs_time`, `pfs_event` | Progression-free survival (months; 1 = progression or death) |

## Requirements

R (>= 4.1) and the following packages:

```r
install.packages(c("survival", "rms", "gtsummary", "survminer",
                   "ggplot2", "dplyr", "scales"))
```

## Usage

From the repository root:

```r
for (f in list.files("R", full.names = TRUE)) source(f)
data <- load_septralu()

table1(data)              # Baseline characteristics by primary tumour site
table2(data)              # Treatment response (RECIST, clinical, biochemical)
table3(data)              # Kaplan-Meier estimates by primary tumour site
figure1(data)             # Overall and progression-free survival curves
figure2(data)             # Forest plot of multivariable Cox models
toxicity_table(data)      # Maximum toxicity grade per patient (CTCAE)
```

Each `table*` function returns a table object or data frame; each `figure*`
function returns a plot. `table2()` also reports the disease control and
objective response rates.

To regenerate everything at once and write the figures to `output/`:

```sh
Rscript reproduce.R
```

## Exploratory analyses

```r
cycles_analysis(data)     # Outcomes by number of R-PRRT cycles (2 vs 4)
cycles_timedep(data)      # Number of cycles as a time-dependent covariate
interval_analysis(data)   # Effect of the I-PRRT to R-PRRT interval (continuous, tertiles, cutpoints)
sequence_analysis(data)   # Outcomes by intercurrent systemic therapy between courses
```

`cycles_analysis()` compares survival and response between patients who received
2 and 4 cycles. Because the number of cycles delivered depends on remaining
progression-free (a patient cannot complete 4 cycles after early progression),
a naive baseline comparison is subject to guarantee-time bias.
`cycles_timedep()` addresses this by modelling the cumulative number of cycles
received as a time-dependent covariate in a counting-process Cox model, and
reports the naive estimate alongside for contrast. `interval_analysis()` models
the interval as a continuous predictor and by tertiles.
