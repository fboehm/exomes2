#!/usr/bin/env Rscript
# Summarize the CpG mutation spectrum from the cpg_spectrum pipeline outputs.
#
# From variant_level.tsv + gene_level.refined.tsv, writes:
#   spectrum_by_context.tsv      folded C-centric spectrum x CpG context (all variants)
#   spectrum_singletons.tsv      same, restricted to singletons (AC==1; de novo proxy)
#   spectrum_by_consequence.tsv  per-variant C>T fraction by VEP consequence x context
#   codon_class.tsv              codon geometry of coding CpGs (codon_class x partner_coding)
#   spectrum.png                 faceted bar plot, all vs singleton proportions
#
# Run via the Snakemake `summarize` rule, or standalone from cpg_spectrum/:
#   Rscript scripts/summarize.R

suppressPackageStartupMessages(library(tidyverse))

# --- resolve I/O: Snakemake `script:` object if present, else defaults --------
if (exists("snakemake")) {
  variant_tsv     <- snakemake@input[["variant"]]
  gene_tsv        <- snakemake@input[["gene"]]
  out_spectrum    <- snakemake@output[["spectrum"]]
  out_singletons  <- snakemake@output[["singletons"]]
  out_consequence <- snakemake@output[["consequence"]]
  out_codon       <- snakemake@output[["codon"]]
  out_plot        <- snakemake@output[["plot"]]
} else {
  variant_tsv     <- "variant_level.tsv"
  gene_tsv        <- "gene_level.refined.tsv"
  out_spectrum    <- "summary/spectrum_by_context.tsv"
  out_singletons  <- "summary/spectrum_singletons.tsv"
  out_consequence <- "summary/spectrum_by_consequence.tsv"
  out_codon       <- "summary/codon_class.tsv"
  out_plot        <- "summary/spectrum.png"
}
dir.create(dirname(out_spectrum), showWarnings = FALSE, recursive = TRUE)

# --- helpers -----------------------------------------------------------------
COMP <- c(A = "T", C = "G", G = "C", T = "A")   # complement, for strand folding

# The pipeline keeps only ref in {C,G}; fold G>N onto the complementary C>N'
# so the spectrum is purely C-centric (C>A, C>G, C>T).
fold_spectrum <- function(df) {
  df |>
    mutate(
      folded  = if_else(ref == "G", paste0("C>", COMP[alt]), paste0("C>", alt)),
      folded  = factor(folded, levels = c("C>A", "C>G", "C>T")),
      context = if_else(is_cpg == 1L, "CpG", "non-CpG")
    )
}

spectrum_table <- function(df) {
  df |>
    count(context, folded, name = "n") |>
    group_by(context) |>
    mutate(prop = n / sum(n)) |>
    ungroup() |>
    arrange(context, folded)
}

# --- variant-level spectrum --------------------------------------------------
v <- read_tsv(
  variant_tsv,
  col_select = c(ref, alt, ac, is_cpg),
  col_types  = cols(ref = "c", alt = "c", ac = "i", is_cpg = "i")
) |>
  fold_spectrum()

spectrum_all <- spectrum_table(v)
write_tsv(spectrum_all, out_spectrum)

spectrum_singletons <- v |>
  filter(ac == 1L) |>
  spectrum_table()
write_tsv(spectrum_singletons, out_singletons)

# --- functional stratification (gene level) ----------------------------------
g <- read_tsv(
  gene_tsv,
  col_select = c(ref, alt, is_cpg, biotype, consequence, codon_class, partner_coding),
  col_types  = cols(.default = "c", is_cpg = "i")
)

# C>T fraction (the methyl-CpG deamination signal) per VEP consequence x context.
by_consequence <- g |>
  mutate(context = if_else(is_cpg == 1L, "CpG", "non-CpG"),
         is_ct   = (ref == "C" & alt == "T") | (ref == "G" & alt == "A")) |>
  group_by(consequence, context) |>
  summarise(n = n(), ct_frac = mean(is_ct), .groups = "drop") |>
  arrange(consequence, context)
write_tsv(by_consequence, out_consequence)

# Codon geometry of coding CpGs: how the partner base / codon split up.
codon_summary <- g |>
  filter(is_cpg == 1L, biotype == "protein_coding") |>
  count(codon_class, partner_coding, name = "n") |>
  arrange(desc(n))
write_tsv(codon_summary, out_codon)

# --- plot: folded spectrum proportions, all vs singletons --------------------
plot_df <- bind_rows(
  spectrum_all        |> mutate(set = "all variants"),
  spectrum_singletons |> mutate(set = "singletons (AC=1)")
)

p <- ggplot(plot_df, aes(folded, prop, fill = context)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~ set) +
  labs(x = "substitution (C-centric)", y = "proportion of variants",
       fill = NULL, title = "CpG mutation spectrum") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(out_plot, p, width = 8, height = 4, dpi = 150)

message("summarize.R: wrote ", out_spectrum, ", ", out_singletons, ", ",
        out_consequence, ", ", out_codon, ", ", out_plot)
