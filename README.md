# SEPTRALU retreatment

R code for the analyses of [177Lu]Lu-DOTATATE retreatment (R-PRRT) in metastatic
neuroendocrine neoplasms within the Spanish SEPTRALU registry.

## Data

The original data are not included in this repository and are available from
the authors upon reasonable request.

## Requirements

R (>= 4.1) and the following packages:

```r
install.packages(c("survival", "rms", "Hmisc", "gtsummary", "survminer",
                   "ggplot2", "dplyr", "scales"))
```

## Usage

From the repository root, with the analysis dataset (`septralu_retreatment.rds`)
in the same folder:

```r
for (f in c("setup.R", "tables.R", "figures.R", "sensitivity.R", "cycles.R")) source(f)
load_dependencies()
data <- load_septralu()

table1(data)              # Baseline characteristics by primary tumour site
table2(data)              # Treatment response
table3(data)              # Kaplan-Meier estimates by primary tumour site
multivariable_cox(data)   # Multivariable Cox models for OS and PFS
figure1(data)             # Overall and progression-free survival
figure2(data)             # Forest plot of the multivariable models
sensitivity_table(data)   # Sensitivity analyses of the multivariable models
cycles_analysis(data)     # Outcomes by number of R-PRRT cycles
cycles_timedep(data)      # Number of cycles as a time-dependent covariate
interval_analysis(data)   # Interval from initial PRRT to R-PRRT
sequence_analysis(data)   # Systemic therapy between the two PRRT courses
```

`Rscript reproduce.R` runs the main analyses and writes the figures to `output/`.
