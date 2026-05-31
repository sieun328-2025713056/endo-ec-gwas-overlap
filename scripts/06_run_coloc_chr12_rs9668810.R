.libPaths(Sys.getenv("R_LIBS_USER"))

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

# Exploratory p-value-based colocalization for chr12 rs9668810-indexed region.
# Inputs must be prepared beforehand as:
#   coloc_input_adeno_allEC_chr12_rs9668810_500kb.full.tsv
#   coloc_input_adeno_EEC_chr12_rs9668810_500kb.full.tsv
# This analysis is exploratory because adenomyosis GWAS beta/SE were not available.

args <- commandArgs(trailingOnly = TRUE)
project <- Sys.getenv("PROJECT", unset = file.path(Sys.getenv("HOME"), "projects", "endo_ecancer"))
coloc_dir <- Sys.getenv("COLOC", unset = file.path(project, "results", "coloc", "rs9668810_region"))
out_dir <- coloc_dir

if (length(args) >= 1 && nzchar(args[1])) coloc_dir <- args[1]
if (length(args) >= 2 && nzchar(args[2])) out_dir <- args[2]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

s_adeno <- 8753 / (8753 + 415718)
s_allEC <- 12906 / (12906 + 108979)
s_EEC   <- 8758 / (8758 + 46126)

run_coloc_pval <- function(file, outcome_prefix, s_outcome, n_outcome, outprefix) {
  infile <- file.path(coloc_dir, file)
  if (!file.exists(infile)) {
    stop("Missing coloc input file: ", infile)
  }

  d <- fread(infile)

  p2_col   <- paste0(outcome_prefix, "_P")
  maf2_col <- paste0(outcome_prefix, "_MAF")

  required_cols <- c("SNP", "AD_P", "AD_MAF", "AD_N", p2_col, maf2_col)
  missing_cols <- setdiff(required_cols, names(d))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", infile, ": ", paste(missing_cols, collapse = ", "))
  }

  d <- d[
    is.finite(AD_P) &
      is.finite(get(p2_col)) &
      is.finite(AD_MAF) &
      is.finite(get(maf2_col)) &
      AD_P > 0 & AD_P <= 1 &
      get(p2_col) > 0 & get(p2_col) <= 1 &
      AD_MAF > 0 & AD_MAF < 0.5 &
      get(maf2_col) > 0 & get(maf2_col) < 0.5
  ]

  d <- unique(d, by = "SNP")

  ds1 <- list(
    pvalues = d$AD_P,
    MAF = d$AD_MAF,
    snp = d$SNP,
    N = median(d$AD_N, na.rm = TRUE),
    type = "cc",
    s = s_adeno
  )

  ds2 <- list(
    pvalues = d[[p2_col]],
    MAF = d[[maf2_col]],
    snp = d$SNP,
    N = n_outcome,
    type = "cc",
    s = s_outcome
  )

  check_dataset(ds1, suffix = "ADENO")
  check_dataset(ds2, suffix = outcome_prefix)

  res <- coloc.abf(ds1, ds2)

  summary_dt <- as.data.table(as.list(res$summary))
  summary_dt[, analysis := outprefix]
  summary_dt[, method := "pvalue_based"]
  setcolorder(summary_dt, c("analysis", "method", setdiff(names(summary_dt), c("analysis", "method"))))
  fwrite(summary_dt, file.path(out_dir, paste0(outprefix, "_coloc_summary.tsv")), sep = "\t")

  res_dt <- as.data.table(res$results)
  fwrite(res_dt, file.path(out_dir, paste0(outprefix, "_coloc_full_results.tsv")), sep = "\t")

  if ("SNP.PP.H4" %in% names(res_dt)) {
    top_dt <- res_dt[order(-SNP.PP.H4)][1:min(30, .N)]
    fwrite(top_dt, file.path(out_dir, paste0(outprefix, "_top_H4_snps.tsv")), sep = "\t")
  }

  cat("\n====================\n")
  cat(outprefix, "\n")
  cat("====================\n")
  print(res$summary)

  if ("SNP.PP.H4" %in% names(res_dt)) {
    cat("\nTop H4 SNPs:\n")
    print(res_dt[order(-SNP.PP.H4)][1:min(10, .N), .(snp, SNP.PP.H4)])
  }

  invisible(res)
}

res_allEC <- run_coloc_pval(
  file = "coloc_input_adeno_allEC_chr12_rs9668810_500kb.full.tsv",
  outcome_prefix = "EC",
  s_outcome = s_allEC,
  n_outcome = 121885,
  outprefix = "adeno_allEC_rs9668810_500kb"
)

res_EEC <- run_coloc_pval(
  file = "coloc_input_adeno_EEC_chr12_rs9668810_500kb.full.tsv",
  outcome_prefix = "EEC",
  s_outcome = s_EEC,
  n_outcome = 54884,
  outprefix = "adeno_EEC_rs9668810_500kb"
)
