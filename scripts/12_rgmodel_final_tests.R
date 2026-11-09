test_rg <- function(file,outcome){
  x <- readRDS(file)

  r <- c(
    ENDO=x$R[1,4],
    ADENO=x$R[2,4],
    NOADENO=x$R[3,4]
  )

  V <- x$V_R[c(3,5,6),c(3,5,6)]

  C <- rbind(
    c(1,-1,0),
    c(1,0,-1)
  )

  d <- C %*% r
  Q <- drop(t(d) %*%
            solve(C %*% V %*% t(C)) %*%
            d)

  omni <- data.frame(
    Outcome=outcome,
    ChiSq=Q,
    df=2,
    P=pchisq(Q,2,lower.tail=FALSE)
  )

  Cp <- rbind(
    "ENDO-ADENO"=c(1,-1,0),
    "ENDO-NOADENO"=c(1,0,-1),
    "ADENO-NOADENO"=c(0,1,-1)
  )

  diff <- drop(Cp %*% r)
  se <- sqrt(diag(Cp %*% V %*% t(Cp)))

  pair <- data.frame(
    Outcome=outcome,
    Contrast=rownames(Cp),
    Difference=diff,
    SE=se,
    Z=diff/se,
    P=2*pnorm(-abs(diff/se))
  )

  list(omni=omni,pair=pair)
}

A <- test_rg("results/rgmodel_allEC.rds","allEC")
E <- test_rg("results/rgmodel_EEC.rds","EEC")

omni <- rbind(A$omni,E$omni)
pair <- rbind(A$pair,E$pair)
pair$Holm_P_6 <- p.adjust(pair$P,method="holm")

print(omni)
print(pair)
