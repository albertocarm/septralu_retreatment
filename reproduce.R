## Reproduce all tables and figures.
## Run from the repository root: Rscript reproduce.R

for (f in c("setup.R", "tables.R", "figures.R", "sensitivity.R", "cycles.R")) source(f)
load_dependencies()

data <- load_septralu()
dir.create("output", showWarnings = FALSE)

t1 <- table1(data)
t2 <- table2(data)
t3 <- table3(data)
mv <- multivariable_cox(data)
sens <- sensitivity_table(data)

print(t3)
print(mv)
print_sensitivity_table(sens)

cyc <- cycles_analysis(data)
cyc_td <- cycles_timedep(data)
int <- interval_analysis(data)

ggplot2::ggsave("output/figure2.pdf", figure2(data, mv), width = 10, height = 3.39,
                device = grDevices::cairo_pdf)
png("output/figure1.png", width = 1200, height = 600, res = 120)
figure1(data)
dev.off()

saveRDS(list(table1 = t1, table2 = t2, table3 = t3, multivariable = mv,
             sensitivity = sens, cycles = cyc, cycles_timedep = cyc_td,
             interval = int),
        "output/results.rds")
message("Done. Figures written to output/.")
