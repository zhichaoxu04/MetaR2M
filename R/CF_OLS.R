#' CF_OLS.R
#'
#' A more efficient two-stage cross-fitted estimation procedure for the R-squared measure in the presence of high-dimensional omics mediators.
#'
#' @param Y Outcome input.
#' @param M Mediator(s) input.
#' @param Covar Covariates input.
#' @param X Exposure input.
#' @param iter.max Maximum number of iterations for (i)SIS and its variants (Default = 3).
#' @param nsis Number of pedictors recuited by (I)SIS (Default = NULL).
#' @param seed Random seed used for sample splitting, and cross-fitting of two subsamples (If the seed = NULL, the data will be evenly split into two halves).
#' @param FDR Indicator to perform FDR control.
#' @param FDRCutoff Cut-off point for FDR control.
#' @param method Methods for the mediators selection: iSIS or HDMT (Default = iSIS).
#'
#' @return A list with the elements:
#' \item{R2M}{Estimated R2-based total mediation effect.}
#'
#' @importFrom stats rnorm residuals lm var qnorm cov
#' @importFrom SIS SIS
#' @importFrom dplyr select
#' @importFrom HDMT null_estimation fdr_est
#'
#' @export
#' @examples
#' p0 <- 20
CF_OLS <- function(Y, M, Covar, X, iter.max=3, nsis=NULL, seed=2024, FDR=FALSE, FDRCutoff=0.2, method="iSIS"){

  # Set Seed for reproducibility
  set.seed(seed)
  startTime <- Sys.time()

  n_Y <- length(Y)
  n_X <- length(X)
  n_M <- length(M)
  if(n_Y!= n_X | n_Y!= n_M | n_X!= n_M){
    stop("Sample size of Y, X, or M does not match.")
  }else{
    n <- n_Y
    d <- ncol(M)
  }

  # Split the sample into two subsamples
  if (is.null(seed)) {
    idx1 <- 1:(nrow(M) * 1/2)
  } else {
    set.seed(seed)
    idx1 <- sample(1:n, ceiling(n/2), replace = FALSE)
  }

  # Check covariates
  if (!is.null(Covar)){
    Covar <- as.matrix(Covar)

    # Function for M ~ COV
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

    # Residuals: M ~ Cov
    # M ~ Cov
    M_res_3 <- apply(M[idx1, ], 2, orth_3)  # M ~ Cov -> M.res - train
    M_res_3 <- apply(M_res_3, 2, scale)

    M_res_4 <- apply(M[-idx1, ], 2, orth_4)  # M ~ Cov -> M.res - test
    M_res_4 <- apply(M_res_4, 2, scale)

    # M.std
    M_res_5 <- apply(M, 2, scale)[idx1, ]
    M_res_6 <- apply(M, 2, scale)[-idx1, ]

    # Residuals: Y ~ Cov
    tdat_3 <- data.frame(y = Y[idx1],
                         cov = Covar[idx1, ])
    f_3 <- stats::lm(y ~ ., data = tdat_3)
    Y_res_3 <- stats::residuals(f_3)     # Y ~ Cov -> Y.res - train
    Y_res_3.std <- as.numeric(scale(Y_res_3, center = T, scale = T))

    tdat_4 <- data.frame(y = Y[-idx1],
                         cov = Covar[-idx1, ])
    f_4 <- stats::lm(y ~ ., data = tdat_4)
    Y_res_4 <- stats::residuals(f_4)    # Y ~ Cov -> Y.res - test
    Y_res_4.std <- as.numeric(scale(Y_res_4, center = T, scale = T))

    # Standardize X
    X_res_3 <- as.numeric(scale(X[idx1]))
    X_res_4 <- as.numeric(scale(X[-idx1]))

    # Residuals: X ~ Cov
    xdat_1 <- data.frame(y = X[idx1],
                         cov = Covar[idx1, ])
    fx_1 <- stats::lm(y ~ ., data = xdat_1)
    X_res_1 <- stats::residuals(fx_1)   # 1/2 n * 1
    X_res_1.std <- as.numeric(scale(X_res_1, center = T, scale = T)) # X ~ Cov -> X.res - train

    xdat_2 <- data.frame(y = X[-idx1],
                         cov = Covar[-idx1, ])
    fx_2 <- stats::lm(y ~ ., data = xdat_2)
    X_res_2 <- stats::residuals(fx_2)   # 1/2 n * 1
    X_res_2.std <- as.numeric(scale(X_res_2, center = T, scale = T)) # X ~ Cov -> X.res - test
  }else{
    # X
    X_res_1 <- X[idx1]
    X_res_1.std <- as.numeric(scale(X_res_1, center = T, scale = T))
    X_res_2 <- X[-idx1]
    X_res_2.std <- as.numeric(scale(X_res_2, center = T, scale = T))
    X_res_3 <- as.numeric(scale(X[idx1]))
    X_res_4 <- as.numeric(scale(X[-idx1]))

    # Y
    Y_res_3 <- Y[idx1]
    Y_res_3.std <- as.numeric(scale(Y_res_3, center = T, scale = T))
    Y_res_4 <- Y[-idx1]
    Y_res_4.std <- as.numeric(scale(Y_res_4, center = T, scale = T))

    # M
    M_res_3 <- M[idx1, ]
    M_res_3 <- apply(M_res_3, 2, scale)

    M_res_4 <- M[-idx1, ]
    M_res_4 <- apply(M_res_4, 2, scale)
  }

  # Regress Y on X
  tdat1 <- data.frame(y = Y_res_3, x = X_res_1)
  f0 <- stats::lm(y~., data=tdat1)
  Y_res_5 <- stats::residuals(f0)

  tdat2 <- data.frame(y = Y_res_4, x = X_res_2)
  f1 <- stats::lm(y~., data=tdat2)
  Y_res_6 <- stats::residuals(f1)

  # Regress M on X
  MX_1 <- function(x){
    data1 <- data.frame(Med = x,
                        envir = X_res_1)
    l <- summary(stats::lm('Med ~.', data = data1))
    res <- stats::residuals(l)
    return(res)
  }

  MX_2 <- function(x){
    data2 <- data.frame(Med = x,
                        envir = X_res_2)
    l <- summary(stats::lm('Med ~.', data = data2))
    res <- stats::residuals(l)
    return(res)
  }

  M_res_1 <- apply(M_res_3, 2, MX_1)
  M_res_1 <- apply(M_res_1, 2, scale)  # M.res ~ X.res - train

  M_res_2 <- apply(M_res_4, 2, MX_2)
  M_res_2 <- apply(M_res_2, 2, scale)  # M.res ~ X.res - test

  Lmc1 <- stats::lm(Y_res_4 ~  X_res_2)
  Allc1 <- Lmc1$coefficients[2]
  AllcPvalue1 <- summary(Lmc1)$coefficients[2, 4]

  Lmc2 <- stats::lm(Y_res_3 ~  X_res_1)
  Allc2 <- Lmc2$coefficients[2]
  AllcPvalue2 <- summary(Lmc2)$coefficients[2, 4]

  cal_alpha_simple <- function(yyy, xxx, type=1){
    data2 <- data.frame(Med = yyy, envir = xxx)
    l <- summary(stats::lm('Med ~ .', data = data2))
    if (type == 1){
      value <- (l$coefficients)['envir','Pr(>|t|)']
    }else{
      value <- as.numeric(l$coefficients['envir','Estimate'])
    }
    return(value)
  }

  cal_beta_simple <- function(yyy, xxx, zzz, type=1){
    data2 <- data.frame(Med = yyy, envir = xxx, outcome = zzz)
    l <- summary(stats::lm('outcome ~ Med + envir', data = data2))
    if (type == 1){
      value <- (l$coefficients)['Med','Pr(>|t|)']
    }else if (type == 2){
      value <- as.numeric(l$coefficients['Med','Estimate'])
    }
    return(value)
  }

  AllAlphaPvalue1 <- apply(M_res_4, 2, function(yyy) cal_alpha_simple(yyy, xxx=X_res_2, type=1))
  AllAlpha1 <- apply(M_res_4, 2, function(yyy) cal_alpha_simple(yyy, xxx=X_res_2, type=2))

  AllAlphaPvalue2 <- apply(M_res_3, 2, function(yyy) cal_alpha_simple(yyy, xxx=X_res_1, type=1))
  AllAlpha2 <- apply(M_res_3, 2, function(yyy) cal_alpha_simple(yyy, xxx=X_res_1, type=2))

  AllBetaPvalue1 <- apply(M_res_4, 2, function(yyy) cal_beta_simple(yyy, xxx=X_res_2, zzz=Y_res_4, type=1))
  AllBeta1 <- apply(M_res_4, 2, function(yyy) cal_beta_simple(yyy, xxx=X_res_2, zzz=Y_res_4, type=2))

  AllBetaPvalue2 <- apply(M_res_3, 2, function(yyy) cal_beta_simple(yyy, xxx=X_res_1, zzz=Y_res_3, type=1))
  AllBeta2 <- apply(M_res_3, 2, function(yyy) cal_beta_simple(yyy, xxx=X_res_1, zzz=Y_res_3, type=2))


  # ------------------- Mediators Selection --------
  if(method == "HDMT"){
    nullEst1 <- HDMT::null_estimation(cbind(AllAlphaPvalue1, AllBetaPvalue1))
    JM_FDR1 <- HDMT::fdr_est(alpha00 = nullEst1$alpha00,
                             alpha01 = nullEst1$alpha01,
                             alpha10 = nullEst1$alpha10,
                             alpha1 = nullEst1$alpha1,
                             alpha2 = nullEst1$alpha2,
                             input_pvalues = cbind(AllAlphaPvalue1, AllBetaPvalue1),
                             exact = 0)
    m1 <- which(JM_FDR1 <= FDRCutoff)
    pabBefore_1 <- length(m1)
    pabAfter_1 <- length(m1)
    nullEst2 <- HDMT::null_estimation(cbind(AllAlphaPvalue2, AllBetaPvalue2))
    JM_FDR2 <- HDMT::fdr_est(alpha00 = nullEst2$alpha00,
                             alpha01 = nullEst2$alpha01,
                             alpha10 = nullEst2$alpha10,
                             alpha1 = nullEst2$alpha1,
                             alpha2 = nullEst2$alpha2,
                             input_pvalues = cbind(AllAlphaPvalue2, AllBetaPvalue2),
                             exact = 0)
    m2 <- which(JM_FDR2 <= FDRCutoff)
    pabBefore_2 <- length(m2)
    pabAfter_2 <- length(m2)
    if(pabAfter_1 == 0 | pabAfter_1 == 0){
      cat("There is no mediators selected.", "\n")
      break
    }

  }else if(method == "iSIS"){
    # iSIS in subsample 2
    model1 <- invisible(SIS::SIS(x = M_res_1, y = Y_res_5,
                                 family = 'gaussian', tune = 'bic',  seed = 1234,
                                 penalty = 'MCP',
                                 nsis = nsis,
                                 iter.max = iter.max))
    pabBefore_1 <- length(model1$ix)
    m1 <- model1$ix
    pabAfter_1 <- length(m1)

    if (length(m1) > 0){
      M1AlphaPvalue <- apply(M_res_4[, m1], 2, function(yyy) cal_alpha_simple(yyy, xxx=X_res_2, type=1)) # M ~ X
      M1Alpha <- apply(M_res_4[, m1], 2, function(yyy) cal_alpha_simple(yyy, xxx=X_res_2, type=2))

      if(FDR == TRUE){
        m1 <- m1[which(M1AlphaPvalue <= FDRCutoff)]
        M1Alpha <- M1Alpha[which(M1AlphaPvalue <= FDRCutoff)]
        M1AlphaPvalue <- M1AlphaPvalue[which(M1AlphaPvalue <= FDRCutoff)]
        pabAfter_1 <- length(m1)
      }
    }else{
      paste0("There is no mediators selected in the 1st half")
    }

    # iSIS in subsample 2
    model2 <- invisible(SIS::SIS(x = M_res_2, y = Y_res_6,
                                 family = 'gaussian', tune = 'bic',  seed = 1234,
                                 penalty = 'MCP',
                                 nsis = nsis,
                                 iter.max = iter.max))
    pabBefore_2 <- length(model2$ix)
    m2 <- model2$ix
    pabAfter_2 <- length(m2)

    if (length(m2) > 0) {
      M2AlphaPvalue <- apply(M_res_3[, m2], 2, function(yyy) cal_alpha_simple(yyy, xxx=X_res_1, type=1)) # M ~ X
      M2Alpha <- apply(M_res_3[, m2], 2, function(yyy) cal_alpha_simple(yyy, xxx=X_res_1, type=2))

      if(FDR == TRUE){
        m2 <- m2[which(M2AlphaPvalue <= FDRCutoff)]
        M2Alpha <- M2Alpha[which(M2AlphaPvalue <= FDRCutoff)]
        M2AlphaPvalue <- M2AlphaPvalue[which(M2AlphaPvalue <= FDRCutoff)]
        pabAfter_2 <- length(m2)
      }
    }else{
      paste0("There is no mediators selected in the 2nd half")
    }
  }

  # Estimation procedure
  ols1.YW <- stats::lm(Y_res_4 ~ cbind(M_res_4, X_res_2.std)[, c(m1, d + 1)])   # Y ~ M + X
  err_yw1 <- ols1.YW$residuals

  M1Beta <- ols1.YW$coefficients[2:(length(m1)+1)]
  M1BetaPvalue <- summary(ols1.YW)$coefficients[, 4][2:(length(m1)+1)]
  M1Gamma <- ols1.YW$coefficients[(length(m1)+2)]
  M1GammaPvalue <- summary(ols1.YW)$coefficients[(length(m1)+2), 4]

  ols1.YZ <- stats::lm(Y_res_4 ~ M_res_4[, m1])                             # Y ~ M
  err_yz1 <- ols1.YZ$residuals

  ols1.YX <- stats::lm(Y_res_4 ~ X_res_2)                                    # Y ~ X
  err_yx1 <- ols1.YX$residuals
  M1c <- ols1.YX$coefficients[2]
  M1cPvalue <- summary(ols1.YX)$coefficients[2, 4]

  # ABproduct_1 <- sum(as.numeric(M1Beta) * as.numeric(M1Alpha))

  err_y1 <- Y_res_4 - mean(Y_res_4)

  v_yx1 <- mean(err_yx1^2)
  v_yw1 <- mean(err_yw1^2)
  v_yz1 <- mean(err_yz1^2)
  v_y1 <- mean(err_y1^2)

  err1 <- cbind(err_yx1^2, err_yz1^2, err_yw1^2, err_y1^2)
  A1 <- stats::cov(err1)

  ols2.YW <- stats::lm(Y_res_3 ~ cbind(M_res_3, X_res_1.std)[, c(m2, d + 1)])
  err_yw2 <- ols2.YW$residuals

  M2Beta <- ols2.YW$coefficients[2:(length(m2)+1)]
  M2BetaPvalue <- summary(ols2.YW)$coefficients[, 4][2:(length(m2)+1)]
  M2Gamma <- ols2.YW$coefficients[(length(m2)+2)]
  M2GammaPvalue <- summary(ols2.YW)$coefficients[(length(m2)+2), 4]

  ols2.YZ <- stats::lm(Y_res_3 ~ M_res_3[, m2])
  err_yz2 <- ols2.YZ$residuals

  ols2.YX <- stats::lm(Y_res_3 ~ X_res_1)
  err_yx2 <- ols2.YX$residuals
  M2c <- ols2.YX$coefficients[2]
  M2cPvalue <- summary(ols2.YX)$coefficients[2, 4]

  # ABproduct_2 <- sum(as.numeric(M2Beta) * as.numeric(M2Alpha))

  err_y2 <- Y_res_3 - mean(Y_res_3)

  v_yx2 <- mean(err_yx2^2)
  v_yw2 <- mean(err_yw2^2)
  v_yz2 <- mean(err_yz2^2)
  v_y2 <- mean(err_y2^2)

  err2 <- cbind(err_yx2^2, err_yz2^2, err_yw2^2, err_y2^2)
  A2 <- stats::cov(err2)

  A <- 0.5 * (A1 + A2)
  v_yw <- 0.5 * (v_yw1 + v_yw2)
  v_yz <- 0.5 * (v_yz1 + v_yz2)
  v_yx <- 0.5 * (v_yx1 + v_yx2)
  v_y <- stats::var(c(Y_res_3, Y_res_4))
  v_y1 <- stats::var(Y_res_4)
  v_y2 <- stats::var(Y_res_3)

  ols.YX <- lm(c(Y_res_3, Y_res_4) ~ c(X_res_3, X_res_4))
  RYX_12 <- summary(ols.YX)$adj.r.squared

  # RYX
  R_YX1 <- 1 - v_yx1 / v_y1
  R_YX2 <- 1 - v_yx2 / v_y2
  R_YX <- 1 - v_yx / v_y
  R_YX_total <- RYX_12

  # RYM
  R_YM1 <- 1 - v_yz1 / v_y1
  R_YM2 <- 1 - v_yz2 / v_y2
  R_YM <- 1 - v_yz / v_y

  # RYXM
  R_YXM1 <- 1 - v_yw1 / v_y1
  R_YXM2 <- 1 - v_yw2 / v_y2
  R_YXM <- 1 - v_yw / v_y

  a <- c(-1/v_y, -1/v_y, 1/v_y, (v_yx + v_yz - v_yw)/v_y^2)
  v <- t(a) %*% A %*% a

  # Compute the R2M
  Rsq.mediated <- 1.0 - (v_yx + v_yz - v_yw) / v_y
  CI_width_asym <- qnorm(0.975) * sqrt(v) / sqrt(n)
  v_asym <- sqrt(v) / sqrt(n)

  # Compute SOS
  SOS <- Rsq.mediated/R_YX
  a_SOS <- c(-(v_yz-v_yw)/(v_y-v_yx)^2, -1/(v_y-v_yx), 1/(v_y-v_yx), (v_yz - v_yw)/(v_y-v_yx)^2)
  v_SOS <- t(a_SOS) %*% A %*% a_SOS
  v_asym_SOS <- sqrt(v_SOS) / sqrt(n)
  SOS_CI_l <- SOS - qnorm(0.975) * sqrt(v_SOS) / sqrt(n)
  SOS_CI_u <- SOS + qnorm(0.975) * sqrt(v_SOS) / sqrt(n)

  endTime <- Sys.time()
  TimeUsed <- difftime(endTime, startTime, units = "mins")
  return(list(output = c(R2M = Rsq.mediated, SE = v_asym,  CI_width_asym = CI_width_asym,
                         CI_low_asym = Rsq.mediated - CI_width_asym, CI_upper_asym = Rsq.mediated + CI_width_asym,
                         pabBefore_1 = round(pabBefore_1, 0), pabAfter_1 = round(pabAfter_1, 0),
                         pabBefore_2 = round(pabBefore_2, 0), pabAfter_2 = round(pabAfter_2, 0),
                         R_YX1=R_YX1, R_YX2=R_YX2, R_YX=R_YX, R_YX_total=R_YX_total,
                         R_YM1=R_YM1, R_YM2=R_YM2, R_YM=R_YM,
                         R_YXM1=R_YXM1, R_YXM2=R_YXM2, R_YXM=R_YXM,
                         SOS=SOS, v_asym_SOS=v_asym_SOS, SOS_CI_l=SOS_CI_l, SOS_CI_u=SOS_CI_u,
                         M1c=as.numeric(M1c), M2c=as.numeric(M2c), M1Gamma=as.numeric(M1Gamma), M2Gamma=as.numeric(M2Gamma),
                         Time = TimeUsed,
                         n = round(n,0), d = round(d,0)),
              select1 = m1,
              select2 = m2,

              AllBeta1 = AllBeta1, AllBetaPvalue1 = AllBetaPvalue1, AllAlpha1=AllAlpha1, AllAlphaPvalue1=AllAlphaPvalue1,
              Allc1 = Allc1, AllcPvalue1=AllcPvalue1,
              AllBeta2 = AllBeta2, AllBetaPvalue2 = AllBetaPvalue2, AllAlpha2=AllAlpha2, AllAlphaPvalue2=AllAlphaPvalue2,
              Allc2 = Allc2, AllcPvalue2=AllcPvalue2))
}

