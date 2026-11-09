# MAF >=0.05 sensitivity: post-hoc reproducibility verification
# Newly prepared verification script; not the original execution record.
# Uses existing four-trait LDSC RDS files. Does not rerun LDSC.
# Requires explicit opt-in. Does not write or overwrite results.

if (Sys.getenv("RUN_MAF05_VERIFICATION") != "1") {
  stop("Set RUN_MAF05_VERIFICATION=1 to enable verification.")
}

lib <- Sys.getenv("GSEM_LIB", "")
if (nzchar(lib)) {
  .libPaths(c(path.expand(lib), .libPaths()))
}

suppressPackageStartupMessages(library(GenomicSEM))

base <- Sys.getenv(
  "MAF_LDSC_DIR",
  file.path(path.expand("~"), "GenomicSEM_PoC_1", "results", "maf05")
)

check_outcome <- function(outcome, cancer_trait) {
  infile <- file.path(base, paste0("ldsc_", outcome, "_4trait.rds"))
  if (!file.exists(infile)) stop("Missing input: ", infile)

  x <- readRDS(infile)

  expected <- c("ENDO", "ADENO", "NOADENO", cancer_trait)
  stopifnot(identical(colnames(x$S), expected))

  y <- rgmodel(x)

  r <- c(
    ENDO = y$R[1,4],
    ADENO = y$R[2,4],
    NOADENO = y$R[3,4]
  )
  V <- y$V_R[c(3,5,6), c(3,5,6), drop=FALSE]

  stopifnot(
    all(is.finite(r)),
    all(is.finite(V)),
    isTRUE(all.equal(V, t(V), tolerance=1e-7))
  )

  C <- rbind(c(1,-1,0), c(1,0,-1))
  d <- drop(C %*% r)
  W <- C %*% V %*% t(C)

  stopifnot(
    all(eigen(W, symmetric=TRUE, only.values=TRUE)$values > 0)
  )

  Q <- drop(t(d) %*% solve(W) %*% d)

  omni <- data.frame(
    Outcome=outcome,
    ChiSq=Q,
    df=2,
    P=pchisq(Q, 2, lower.tail=FALSE)
  )

  Cp <- rbind(
    "ENDO-ADENO"=c(1,-1,0),
    "ENDO-NOADENO"=c(1,0,-1),
    "ADENO-NOADENO"=c(0,1,-1)
  )

  delta <- drop(Cp %*% r)
  se <- sqrt(diag(Cp %*% V %*% t(Cp)))

  pair <- data.frame(
    Outcome=outcome,
    Contrast=rownames(Cp),
    Difference=delta,
    SE=se,
    P=2*pnorm(-abs(delta/se))
  )

  list(omni=omni, pair=pair)
}

a <- check_outcome("allEC", "ALLEC")
e <- check_outcome("EEC", "EEC")

omni <- rbind(a$omni, e$omni)
pair <- rbind(a$pair, e$pair)
pair$Holm_P_6 <- p.adjust(pair$P, method="holm")

cat("\n=== MAF sensitivity omnibus ===\n")
print(omni)

cat("\n=== MAF sensitivity pairwise ===\n")
print(pair)

# Previously reported values: comparison only, not inputs to computation.
reported <- c(
  allEC_omni=0.02458479,
  EEC_omni=0.01214624,
  allEC_ENDO_ADENO_Holm=0.08779310,
  EEC_ENDO_ADENO_Holm=0.05680075
)

computed <- c(
  allEC_omni=omni$P[omni$Outcome=="allEC"],
  EEC_omni=omni$P[omni$Outcome=="EEC"],
  allEC_ENDO_ADENO_Holm=pair$Holm_P_6[
    pair$Outcome=="allEC" & pair$Contrast=="ENDO-ADENO"
  ],
  EEC_ENDO_ADENO_Holm=pair$Holm_P_6[
    pair$Outcome=="EEC" & pair$Contrast=="ENDO-ADENO"
  ]
)

cat("\n=== Comparison with previously reported values ===\n")
print(data.frame(
  Metric=names(reported),
  Reported=unname(reported),
  Computed=unname(computed),
  Absolute_difference=abs(computed-reported)
))

if (any(abs(computed-reported) > 1e-5)) {
  stop("Reported and computed results differ; investigate.")
}

cat("\nVERIFICATION PASSED\n")
