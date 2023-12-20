#' CF_OLS.R
#'
#' Calculate the test statistic and p-value for interval-censored competing risks SKAT.
#'
#' @param leftDmat n*(p+nknots+2) design matrix for left end of
#'
#' @return A list with the elements:
#'
#' @importFrom stats rnorm
#'
#' @export
#' @examples
#' p0 <- 20
CF_OLS <- function(Y, M, Covar, X, d, n, iter.max = 3, nsis = NULL, first.half = TRUE, seed = 2022, idx1 = NULL){
  # ----- Load the packages
  library(GMMAT)
  library(SIS)
  library(tidyverse)

  # Set Seed for producibility
  set.seed(seed)
  startTime <- Sys.time()

  # ------ iSIS for mediators selection ------
  # Subset 1 Estimation
  if(!is.null(idx1)){
    idx1 <- idx1
  }else{
    if (first.half == TRUE) {
      idx1 <- 1:(nrow(M) * 1/2)
    } else {
      set.seed(seed)
      idx1 <- sample(1:n, ceiling(n/2), replace = FALSE)
    }

  }


  Covar <- as.matrix(Covar)
  # if(is.null(nsis)){nsis <- length(Y)/log(length(Y))}

  # Residuas M ~ Cov + X.std
  orth_1 <- function(x){

    data1 <- data.frame(Med = x,
                        envir = scale(X[idx1]),
                        cov = Covar[idx1, ])

    l <- summary(stats::lm('Med ~.', data = data1))
    res <- stats::residuals(l)
    return(res)
  }

  orth_2 <- function(x){
    data2 <- data.frame(Med = x,
                        envir = scale(X[-idx1]),
                        cov = Covar[-idx1, ])
    l <- summary(stats::lm('Med ~.', data = data2))
    res <- stats::residuals(l)
    return(res)
  }

  orth_3 <- function(x){

    data3 <- data.frame(Med = x,
                        # envir = scale(X[idx1]),
                        cov = Covar[idx1, ])

    l <- summary(stats::lm('Med ~.', data = data3))
    res <- stats::residuals(l)
    return(res)
  }

  orth_4 <- function(x){
    data4 <- data.frame(Med = x,
                        # envir = scale(X[-idx1]),
                        cov = Covar[-idx1, ])
    l <- summary(stats::lm('Med ~.', data = data4))
    res <- stats::residuals(l)
    return(res)
  }


  # 1st residuals: M ~ Cov + X.std
  M_res_1 <- apply(M[idx1, ], 2, orth_1)  # 1/2 n * d
  M_res_1 <- apply(M_res_1, 2, scale)

  M_res_2 <- apply(M[-idx1, ], 2, orth_2)  # 1/2 n * d
  M_res_2 <- apply(M_res_2, 2, scale)

  # M ~ Cov
  M_res_3 <- apply(M[idx1, ], 2, orth_3)  # 1/2 n * d
  M_res_3 <- apply(M_res_3, 2, scale)

  M_res_4 <- apply(M[-idx1, ], 2, orth_4)  # 1/2 n * d
  M_res_4 <- apply(M_res_4, 2, scale)

  # M.std
  M_res_5 <- apply(M, 2, scale)[idx1, ]
  M_res_6 <- apply(M, 2, scale)[-idx1, ]

  # 2nd residuals: Y ~ Cov + X.std
  tdat_1 <- data.frame(y = Y[idx1],
                       envir = scale(X[idx1]),
                       cov = Covar[idx1, ])
  f_1 <- stats::lm(y ~ ., data = tdat_1)
  Y_res_1 <- stats::residuals(f_1)

  tdat_2 <- data.frame(y = Y[-idx1],
                       envir = scale(X[-idx1]),
                       cov = Covar[-idx1, ])
  f_2 <- stats::lm(y ~ ., data = tdat_2)
  Y_res_2 <- stats::residuals(f_2)

  # 3rd residuals: Y ~ Cov
  tdat_3 <- data.frame(y = Y[idx1],
                       cov = Covar[idx1, ])
  f_3 <- stats::lm(y ~ ., data = tdat_3)
  Y_res_3 <- stats::residuals(f_3)

  tdat_4 <- data.frame(y = Y[-idx1],
                       cov = Covar[-idx1, ])
  f_4 <- stats::lm(y ~ ., data = tdat_4)
  Y_res_4 <- stats::residuals(f_4)

  # 4th x.std
  X_res_3 <- as.numeric(scale(X[idx1]))
  X_res_4 <- as.numeric(scale(X[-idx1]))


  # 4th residuals: X ~ Cov
  # xdat_1 <- data.frame(y = X[idx1],
  #                    cov = Covar[idx1, ])
  # fx_1 <- stats::lm(y ~ ., data = xdat_1)
  # X_res_1 <- stats::residuals(fx_1)   # 1/2 n * 1
  # X_res_1 <- as.numeric(scale(X_res_1, center = T, scale = T))
  #
  # xdat_2 <- data.frame(y = X[-idx1],
  #                    cov = Covar[-idx1, ])
  # fx_2 <- stats::lm(y ~ ., data = xdat_2)
  # X_res_2 <- stats::residuals(fx_2)   # 1/2 n * 1
  # X_res_2 <- as.numeric(scale(X_res_2, center = T, scale = T))

  # ------------------- Mediators Selection --------

  # ------ iSIS
  model1 <- invisible(SIS::SIS(x = M_res_5, y = Y_res_1,
                               family = 'gaussian', tune = 'bic',  seed = 1234,
                               penalty = 'MCP',
                               nsis = nsis,
                               iter.max = iter.max))
  pab_1 <- length(model1$ix)
  m1 <- model1$ix

  # ------ Estimation Procedure ------
  if (length(m1) > 0){
    ols1.YW <- lm(Y_res_4 ~ cbind(M_res_6, X_res_4)[, c(m1, d + 1)])   # Y ~ M + X
    err_yw1 <- ols1.YW$residuals

    ols1.YZ <- lm(Y_res_4 ~ M_res_6[, m1])                             # Y ~ M
    err_yz1 <- ols1.YZ$residuals

    ols1.YX <- lm(Y_res_4 ~ X_res_4)                                    # Y ~ X
    err_yx1 <- ols1.YX$residuals

    RYX_1 <- summary(ols1.YX)$adj.r.squared
    V.YW1 <- sum(ols1.YW$residuals^2)/ols1.YW$df.residual
    V.YZ1 <- sum(ols1.YZ$residuals^2)/ols1.YZ$df.residual
    V.YX1 <- sum(ols1.YX$residuals^2)/ols1.YX$df.residual

    err_y1 <- Y_res_4 - mean(Y_res_4)

    v_yx1 <- mean(err_yx1^2)
    v_yw1 <- mean(err_yw1^2)
    v_yz1 <- mean(err_yz1^2)
    v_y1 <- mean(err_y1^2)

    err1 <- cbind(err_yx1^2, err_yz1^2, err_yw1^2, err_y1^2)
    A1 <- cov(err1)


  }else{
    paste0("There is no mediators selected in the 1st half")
  }

  # ------ Subset 2 Estimation ------

  # ------ iSIS
  model2 <- invisible(SIS::SIS(x = M_res_6, y = Y_res_2,
                               family = 'gaussian', tune = 'bic',  seed = 1234,
                               penalty = 'MCP',
                               nsis = nsis,
                               iter.max = iter.max))
  pab_2 <- length(model2$ix)
  m2 <- model2$ix


  # -----Estimation Process
  if (length(m2) > 0) {
    ols2.YW <- lm(Y_res_3 ~ cbind(M_res_5, X_res_3)[, c(m2, d + 1)])
    err_yw2 <- ols2.YW$residuals

    ols2.YZ <- lm(Y_res_3 ~ M_res_5[, m2])
    err_yz2 <- ols2.YZ$residuals

    ols2.YX <- lm(Y_res_3 ~ X_res_3)
    err_yx2 <- ols2.YX$residuals

    RYX_2 <- summary(ols2.YX)$adj.r.squared
    V.YW2 <- sum(ols2.YW$residuals^2)/ols2.YW$df.residual
    V.YZ2 <- sum(ols2.YZ$residuals^2)/ols2.YZ$df.residual
    V.YX2 <- sum(ols2.YX$residuals^2)/ols2.YX$df.residual

    err_y2 <- Y_res_3 - mean(Y_res_3)

    v_yx2 <- mean(err_yx2^2)
    v_yw2 <- mean(err_yw2^2)
    v_yz2 <- mean(err_yz2^2)
    v_y2 <- mean(err_y2^2)

    err2 <- cbind(err_yx2^2, err_yz2^2, err_yw2^2, err_y2^2)
    A2 <- cov(err2)

  }else{
    paste0("There is no mediators selected in the 2nd half")
  }

  V.YW <- 0.5 * (V.YW1 + V.YW2)
  V.YZ <- 0.5 * (V.YZ1 + V.YZ2)

  A <- 0.5 * (A1 + A2)
  v_yw <- 0.5 * (v_yw1 + v_yw2)
  v_yz <- 0.5 * (v_yz1 + v_yz2)
  v_yx <- 0.5 * (v_yx1 + v_yx2)
  v_y <- var(c(Y_res_3, Y_res_4))

  a <- c(-1/v_y, -1/v_y, 1/v_y, (v_yx + v_yz - v_yw)/v_y^2)
  v <- t(a) %*% A %*% a

  V.YW <- 0.5 * (V.YW1 + V.YW2)
  V.YZ <- 0.5 * (V.YZ1 + V.YZ2)

  # Compute the R2
  Rsq.mediated <- 1.0 - (v_yx + v_yz - v_yw) / v_y
  CI_width_asym <- qnorm(0.975) * sqrt(v) / sqrt(n)
  v_asym <- sqrt(v) / sqrt(n)


  # V.YX <- 0.5 * (V.YX1 + V.YX2) # 1
  ols.YX <- lm(c(Y_res_3, Y_res_4) ~ c(X_res_3, X_res_4))
  V.YX <- sum(ols.YX$residuals^2)/ols.YX$df.residual
  RYX_12 <- summary(ols.YX)$adj.r.squared
  V.Y  <- var(c(Y_res_3, Y_res_4))

  endTime <- Sys.time()
  TimeUsed <- difftime(endTime, startTime, units = "mins")


  return(list(output = c(Rsq.med = Rsq.mediated,  CI_width_asym = CI_width_asym,
                         CI_low_asym = Rsq.mediated - CI_width_asym, CI_upper_asym = Rsq.mediated + CI_width_asym,
                         pab1 = round(pab_1, 0), pab2 = round(pab_2, 0),
                         VYXM = V.YW, VYXM_1 = V.YW1, VYXM_2 = V.YW2, RYXM = 1 - V.YW/V.Y, RYXM_1 = 1 - V.YW1/var(Y_res_4), RYXM_2 = 1 - V.YW2/var(Y_res_3),
                         VYM = V.YZ, VYM_1 = V.YZ1, VYM_2 = V.YZ2, RYM = 1 - V.YZ/V.Y, RYM_1 = 1 - V.YZ1/var(Y_res_4), RYM_2 = 1 - V.YZ2/var(Y_res_3),
                         RYX = 1 - V.YX/V.Y, RYX_1 = RYX_1, RYX_2 = RYX_2, RYX_c = RYX_12,
                         RYX_1_r = 1 - V.YX1/var(Y_res_4), RYX_2_r = 1 - V.YX2/var(Y_res_3),
                         VY = V.Y, VY_1 = var(Y_res_4), VY_2 = var(Y_res_3),
                         v_asym = v_asym, Time = TimeUsed,
                         n = round(n,0), d = round(d,0)),
              select1 = m1,
              select2 = m2))
}

