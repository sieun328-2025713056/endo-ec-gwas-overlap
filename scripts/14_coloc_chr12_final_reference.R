.libPaths(c(file.path(path.expand("~"), "Rlibs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

project <- Sys.getenv(
  "PROJECT",
  unset = file.path(path.expand("~"), "projects", "endo_ecancer")
)

loc <- Sys.getenv(
  "COLOC_INPUT_DIR",
  unset = file.path(project, "results", "coloc", "rs9668810_region")
)

if (!dir.exists(loc)) {
  stop("Colocalization input directory does not exist: ", loc)
}

for (o in c("allEC", "EEC")) {

  tag <- if (o == "allEC") "EC" else "EEC"

  d <- fread(file.path(
    loc,
    paste0("coloc_input_adeno_", o,
           "_chr12_rs9668810_500kb.full.tsv")
  ))

  pc <- paste0(tag, "_P")
  mc <- paste0(tag, "_MAF")
  bc <- paste0(tag, "_BETA_ALIGNED")
  sc <- paste0(tag, "_SE")

  # Keep the original SNP eligibility rules
  d <- d[
    is.finite(AD_P) &
    is.finite(get(pc)) &
    is.finite(AD_MAF) &
    is.finite(get(mc)) &
    AD_P > 0 & AD_P <= 1 &
    get(pc) > 0 & get(pc) <= 1 &
    AD_MAF > 0 & AD_MAF < 0.5 &
    get(mc) > 0 & get(mc) < 0.5
  ]

  d <- unique(d, by = "SNP")

  stopifnot(
    all(is.finite(d[[bc]])),
    all(is.finite(d[[sc]]) & d[[sc]] > 0)
  )

  # Effective-N assumption for ADENO
  ad <- list(
    pvalues = d$AD_P,
    MAF = d$AD_MAF,
    snp = d$SNP,
    N = median(d$AD_N),
    type = "cc",
    s = 0.5
  )

  # ECAC beta / variance
  ec <- list(
    beta = d[[bc]],
    varbeta = d[[sc]]^2,
    snp = d$SNP,
    N = if (o == "allEC") 121885 else 54884,
    type = "cc",
    s = if (o == "allEC") {
      12906 / 121885
    } else {
      8758 / 54884
    }
  )

  fit <- coloc.abf(
    ad, ec,
    p1 = 1e-4,
    p2 = 1e-4,
    p12 = 1e-5
  )

  # Optional SNP-level posterior export.
  # Requires explicit opt-in and never overwrites existing files.
  if (identical(Sys.getenv("EXPORT_COLOC_POSTERIOR", "0"), "1")) {
    outdir <- Sys.getenv("COLOC_POSTERIOR_DIR", "")

    if (!nzchar(outdir) || !dir.exists(outdir)) {
      stop("Set COLOC_POSTERIOR_DIR to an existing directory.")
    }

    dest <- file.path(
      outdir,
      paste0("adeno_", o, "_final_snp_posterior.tsv")
    )

    if (file.exists(dest)) {
      stop("Refusing to overwrite: ", dest)
    }

    fwrite(as.data.table(fit$results), dest, sep = "\t")
    cat("Saved posterior:", dest, "\n")
  }

  top <- fit$results[
    which.max(fit$results$SNP.PP.H4),
  ]

  cat("\n==========", o, "==========\n")
  print(fit$summary)
  cat("ADENO_N:", ad$N, "\n")
  cat("ADENO_s:", ad$s, "\n")
  cat("Top SNP:", top$snp, "\n")
  cat("Top SNP.PP.H4:", top$SNP.PP.H4, "\n")
}
