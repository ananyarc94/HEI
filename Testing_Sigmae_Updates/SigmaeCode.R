# Mathew McLean
# 2014-10-20
# code snippet of parts of sampler relevant for Sigmae update

sum.weight2recall <- sum(weight*recall)

# Get indices for unconstrained parts of V

od.ind <- row(start[[1]]$V) > col(start[[1]]$V)

# for i=4,6,8,10,12 update Vij for 1 <= j <= i-2 # 2+4+6+8+10 = 30 values
# for i=13,...,19, update Vij for 1 <= j <= i-1  # 12 + ... + 18 = 105  values
# update Vii for all i in IND  # length(IND)=13 values
k <- length(IND)
Vii.ind <- matrix(IND, k, 2)
Vij.ind <- matrix(nrow = 30 + 105, ncol = 2)
k <- 1
for (i in seq(4, 12, by = 2)){
  for (j in seq_len(i-2)){
    Vij.ind[k,] <- c(i, j)
    k <- k + 1
  }
}
for (i in 13:19){
  for (j in seq_len(i - 1)){
    Vij.ind[k,] <- c(i, j)
    k <- k + 1
  }
}


What1 <- Xtildei1%*%beta + Utildei
# What2[-one2n2,] is 0
What2[one2n2, ] <- Xtildei2[one2n2, ]%*%beta + Utildei[one2n2, ]
Wresid1  <- Wtildei[, , 1] - What1
Wresid2 <- Wtildei[one2n2, , 2] - What2[one2n2, ]

tempMat <- crossprod(weightMat*Wresid1, Wresid1) +
	                         crossprod(weightMat[one2n2, ]*Wresid2, Wresid2)

V <- UpdateV(sum.weight2recall, r, theta, V, tempMat, Vii.ind, acceptanceEnv,
	     hyperparameters$diag.V, TRUE)
V <- UpdateV(sum.weight2recall, r, theta, V, tempMat, Vij.ind, acceptanceEnv,
	     hyperparameters$offdiag.V, FALSE)
V <- MakePosDef(V, r, theta)

Sigmae <- tcrossprod(V)
iSigmae <- solve(Sigmae)

#' Update one value on diagonal of matrix V
#'
#' Metropolis update for diagonals or off-diagonals of matrix V which satistifies
#' \eqn{\Sigma_{\epsilon}=V^TV} where
#' \eqn{\Sigma_{\epsilon}} is the measurement error variance.
#' Form of update given on page 1481 of Zhang et al (2011).
#' @param sum_weightrecall sum of the product of survey weights and # of recalls
#' @param r current value of r
#' @param theta current value of theta
#' @param V matrix V
#' @param tempMat matrix used by CalcGArg
#' @param ind index of which elements to update on diagonal of V
#' @param aEnv environment containing total number of accepted proposals for R, theta
#'   and V.
#' @param V.hyp length-two vector giving endpoints for random uniform prior for
#' unconstrained elements on the diagonal or off-diagonal of V (depending on which is
#' to be updated)
#' @param diag.update logical; are diagonal or off-diagonal elements of V to be updated?
#' @return New value for entry (i,i) entry of V determine by Metropolis step
#' @note Internal function only called by `ZhangSampler`.  Also updates acceptance
#' rates using `<<-` operator.
#' @keywords internal
UpdateV <- function(sum.weight2recall, r, theta, V, tempMat, inds, aEnv, V.hyp,
                    diag.update){
#  ind <- c(seq(2, 12, by =2), 13:19)
  len <- nrow(inds)

  # diagV <- diag(V)
  # diag.Vcurr <- diagV[ind]
  # diag.Vcand <- diag.Vcurr + 0.4*(runif(DIM) - 0.5)
  Vcand <- runif(len, V[inds] - .2, V[inds] + .2)
  GofSigmaecurr <- CalcGArg(r, theta, V, tempMat)
  for (i in which(Vcand < V.hyp[2] & Vcand > V.hyp[1])){
    ind <- inds[i, , drop = FALSE]
    oldVij <- V[ind]
    V[ind] <- Vcand[i]
    GofSigmaecand <- tryCatch(CalcGArg(r, theta, V, tempMat), error = function(e) -Inf)
    gg <- GofSigmaecand - GofSigmaecurr
    gg <- if (diag.update)
            exp(gg - 0.5*sum.weight2recall * log(Vcand[i]^2) +
              0.5*sum.weight2recall * log(V[ind]^2))
          else
            exp(gg)
    gg <- pmin(1, gg)
    accepted <- runif(1) <= gg
    if (!is.na(accepted)  && accepted){
      aEnv$V[ind] <- aEnv$V[ind] + 1
      GofSigmaecurr <- GofSigmaecand
    }else
      V[ind] <- oldVij
  }
  V
}

MakePosDef <- function(V, r, theta){
  ## V[1, 1] <- 1 # not in MATLAB code
  ## V[2, 1] <- 0 # ""
  V[3, 1] <- r[1] * sin(theta[1])
  V[3, 2] <- r[1] * cos(theta[1])

  # p. 1494 of Zhang et al.: k = q from paper, j = p
  for (k in 2:5){
    V[2*k+1, 1] <- r[k] * sin(theta[(k-1)^2+1])
    for (j in 2:(2*k-1)){
       V[2*k+1, j] <- r[k] * sin(theta[(k-1)^2+j])
       for (jj in seq_len(j-1))
         V[2*k+1, j] <- V[2*k+1, j] * cos(theta[(k-1)^2+jj])
    }

    V[2*k+1, 2*k] <- r[k]
    for (jj in seq_len(2*k-1))
      V[2*k+1, 2*k] <- V[2*k+1, 2*k]*cos(theta[(k-1)^2+jj])
  }
  tind <- seq(3, 11, by=2)
  diag(V)[tind] <- sqrt(1-r^2)

  for (k in tind)
    V[k+1, k] <- -sum(V[k, seq_len(k-1)]*V[k+1,seq_len(k-1)])/ V[k, k]
  V
}
