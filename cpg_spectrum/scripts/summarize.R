#!/usr/bin/env Rscript
# Summarize the CpG mutation spectrum from the cpg_spectrum pipeline outputs.
#
# The inputs are genome-wide and large, so the heavy reads are CHUNKED
# (readr::read_tsv_chunked) and aggregated incrementally: peak memory is bounded
# by CHUNK rows, not by total file size. Each chunk is reduced to a small
# per-category count table and those partials are summed at the end.
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

CHUNK <- 1e6   # rows per chunk; raise if you have RAM to spare, lower if tight

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

COMP <- c(A = "T", C = "G", G = "C", T = "A")   # complement, for strand folding

# Fold G>N onto the complementary C>N' so the spectrum is C-centric
# (the pipeline keeps only ref in {C,G}).
fold_cols <- function(df) {
  df |>
    mutate(
      folded  = if_else(ref == "G", paste0("C>", COMP[alt]), paste0("C>", alt)),
      context = if_else(is_cpg == 1L, "CpG", "non-CpG")
    )
}

# --- variant-level spectrum (chunked) ----------------------------------------
v_types <- cols(.default = col_character(), ac = col_integer(), is_cpg = col_integer())

# Per chunk: counts of folded x context, both for all variants and singletons.
v_chunk <- function(chunk, pos) {
  chunk <- fold_cols(chunk)
  bind_rows(
    chunk |> count(context, folded, name = "n") |> mutate(set = "all variants"),
    chunk |> filter(ac == 1L) |> count(context, folded, name = "n") |>
      mutate(set = "singletons (AC=1)")
  )
}

v_partials <- read_tsv_chunked(
  variant_tsv, DataFrameCallback$new(v_chunk),
  chunk_size = CHUNK, col_types = v_types
)

spectrum_long <- v_partials |>
  group_by(set, context, folded) |>
  summarise(n = sum(n), .groups = "drop") |>
  group_by(set, context) |>
  mutate(prop = n / sum(n)) |>                     # proportion within each context
  ungroup() |>
  mutate(folded = factor(folded, levels = c("C>A", "C>G", "C>T"))) |>
  arrange(set, context, folded)

write_tsv(spectrum_long |> filter(set == "all variants") |> select(-set), out_spectrum)
write_tsv(spectrum_long |> filter(set == "singletons (AC=1)") |> select(-set), out_singletons)

# --- functional stratification (chunked) -------------------------------------
g_types <- cols(.default = col_character(), is_cpg = col_integer())

# Per chunk: (a) consequence x context counts + C>T counts, and (b) coding-CpG
# codon geometry. Tagged with `metric` so one pass yields both summaries.
g_chunk <- function(chunk, pos) {
  chunk <- chunk |>
    mutate(context = if_else(is_cpg == 1L, "CpG", "non-CpG"),
           is_ct   = (ref == "C" & alt == "T") | (ref == "G" & alt == "A"))
  cons <- chunk |>
    group_by(consequence, context) |>
    summarise(n = n(), n_ct = sum(is_ct), .groups = "drop") |>
    mutate(metric = "consequence")
  codon <- chunk |>
    filter(is_cpg == 1L, biotype == "protein_coding") |>
    count(codon_class, partner_coding, name = "n") |>
    mutate(metric = "codon")
  bind_rows(cons, codon)
}

g_partials <- read_tsv_chunked(
  gene_tsv, DataFrameCallback$new(g_chunk),
  chunk_size = CHUNK, col_types = g_types
)

by_consequence <- g_partials |>
  filter(metric == "consequence") |>
  group_by(consequence, context) |>
  summarise(n = sum(n), n_ct = sum(n_ct), .groups = "drop") |>
  mutate(ct_frac = n_ct / n) |>
  select(consequence, context, n, ct_frac) |>
  arrange(consequence, context)
write_tsv(by_consequence, out_consequence)

codon_summary <- g_partials |>
  filter(metric == "codon") |>
  group_by(codon_class, partner_coding) |>
  summarise(n = sum(n), .groups = "drop") |>
  arrange(desc(n))
write_tsv(codon_summary, out_codon)

# --- plot: folded spectrum proportions, all vs singletons --------------------
p <- ggplot(spectrum_long, aes(folded, prop, fill = context)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~ set) +
  labs(x = "substitution (C-centric)", y = "proportion of variants",
       fill = NULL, title = "CpG mutation spectrum") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(out_plot, p, width = 8, height = 4, dpi = 150)

message("summarize.R: wrote ", out_spectrum, ", ", out_singletons, ", ",
        out_consequence, ", ", out_codon, ", ", out_plot)
