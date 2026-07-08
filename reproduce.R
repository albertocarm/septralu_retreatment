## Reproduce all tables and figures.
## Run from the repository root: Rscript reproduce.R

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
load_dependencies()

data <- load_septralu()
dir.create("output", showWarnings = FALSE)

t1 <- table1(data)
t2 <- table2(data)
t3 <- table3(data)
tox <- toxicity_table(data)

print(t3)
print(tox)

cyc <- cycles_analysis(data)
int <- interval_analysis(data)

ggplot2::ggsave("output/figure2.png", figure2(data), width = 10, height = 5, dpi = 300)
png("output/figure1.png", width = 1200, height = 600, res = 120)
figure1(data)
dev.off()

saveRDS(list(table1 = t1, table2 = t2, table3 = t3, toxicity = tox,
             cycles = cyc, interval = int), "output/results.rds")
message("Done. Figures written to output/.")
