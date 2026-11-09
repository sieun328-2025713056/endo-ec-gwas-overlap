.libPaths(c("~/R/gsem019", .libPaths()))
library(GenomicSEM)

run <- function(infile, outfile, outcome, blocks){

  x <- readRDS(infile)
  y <- rgmodel(x)
  # Read-only verification: preserve existing rgmodel RDS files.

  r <- c(y$R[1,4], y$R[2,4], y$R[3,4])
  V <- y$V_R[c(3,5,6), c(3,5,6)]

  C <- rbind(
    c(1,-1,0),
    c(1,0,-1)
  )

  d <- C %*% r
  Q <- drop(
    t(d) %*%
    solve(C %*% V %*% t(C)) %*%
    d
  )

  data.frame(
    Outcome=outcome,
    Blocks=blocks,
    ChiSq=Q,
    df=2,
    P=pchisq(Q,2,lower.tail=FALSE)
  )
}

res <- rbind(
  run("results/ldsc_allEC_b100.rds",
      "results/rgmodel_allEC_b100.rds","allEC",100),

  run("results/ldsc_allEC_b400.rds",
      "results/rgmodel_allEC_b400.rds","allEC",400),

  run("results/ldsc_EEC_b100.rds",
      "results/rgmodel_EEC_b100.rds","EEC",100),

  run("results/ldsc_EEC_b400.rds",
      "results/rgmodel_EEC_b400.rds","EEC",400)
)

print(res)
