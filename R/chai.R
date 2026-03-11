#' Conditional Hypothesis testing using Auxiliary Information (chai) main function
#' @param z A numeric vector that saving the z-statistics.
#' @param X A numeric matrix/vector (covariates/side-information) that corresponding to z.
#' @param K_vec An integer value or a range specifying the numbers of mixture components. Suggest use 2 or above. If there is a vector/range, the model will automatically select the "best" based on BIC. The default is K_vec = 2:6.
#' @param R An integer value indicating the total number of samples the model generate to estimate the \eqn{\pi_0(x)}. The default is R = 100. More details please see the original paper.
#' @return A list with the following components:
#' \item{z}{The input z-statistics.}
#' \item{X}{The input X (side information).}
#' \item{K}{The optimal number of mixture components.}
#' \item{R}{The input R value.}
#' \item{clFDR}{The conditional local FDR  of every hypotheses.}
#' \item{pi0}{The estimated \eqn{\pi_0(x)} of every hypotheses.}
#' \item{post_w}{The posterior weight \eqn{w_{ik}} of every hypotheses (\eqn{i}) belong to which component (\eqn{k}).}
#' \item{ord}{The order of Hypothesis indices by increasing conditional local FDR (from smallest to largest).}
#' \item{clFDR_sorted}{The sorted conditional local FDR (from smallest to largest).
#' @export


chai <- function(z, X, K_vec = 2:6, R = 100) {
  # require(mclust); require(locfdr); require(admix); require(mvtnorm)

  df <- data.frame(as.data.frame(X))
  names(df) <- paste0("x", seq_len(ncol(df)))   # To avoid there is a 'z' column name inside the X
  df$z <- z
  xcols <- setdiff(names(df), "z")

  fit <- mclust::Mclust(df, G = K_vec)   # modelNames = "VVV"
  zMat <- fit$z
  n <- nrow(df)

  lfdr_naive_all <- numeric(n)
  minRatioAll <- numeric(n)
  ratioZ_givenX_all <- numeric(n)
  post_w <- matrix(nrow = n, ncol = fit$G)

  for (i in 1:n) {
    np <- naiveRemoveOneObs(i, df, fit, zMat)
    xVec <- as.numeric(df[i, xcols])
    zVal <- df$z[i]
    cp <- conditionalParamsForX_custom(xVec, np$pi, np$mu, np$Sigma)
    post_w[i,] <- cp$post_weights

    set.seed(123)
    rMix1 <- rGaussianMix(n = R, cp$post_weights, cp$cond_means, sqrt(cp$cond_vars))

    admixMod <- admix::admix_model(knownComp_dist = "norm",
                                   knownComp_param = c("mean" = 0, "sd" = 1))

    result <- admix::admix_estim(samples = list(rMix1),
                                 admixMod = list(admixMod),
                                 est_method = 'PS')

    minRatioAll[i] <- 1 - get_mixing_weights(result)  # pi0

    ratioZ_givenX_all[i] <- ratio_z_given_x_custom(zVal, cp$post_weights, cp$cond_means, cp$cond_vars)

    lfdr_naive_all[i] <- minRatioAll[i] / ratioZ_givenX_all[i]
  }

  ord <- order(lfdr_naive_all)
  lFDR_sorted <- lfdr_naive_all[ord]
  avgFDR <- cumsum(lFDR_sorted) / seq_along(lFDR_sorted)

  return(list(
    z = z,
    X = X,
    K = fit$G,
    R = R,
    clFDR = lfdr_naive_all,
    pi0 = minRatioAll,
    post_w = post_w,
    ord = ord,
    clFDR_sorted = lFDR_sorted,
    avgFDR = avgFDR
  ))
}
