#' R2_Meta.R
#'
#' A novel meta-analysis framework of the R-squared (R2)-based mediation effect estimation for high-dimensional omics mediators.
#'
#' @param Effects Estimated R2M for each study.
#' @param Study An optional vector with study labels (Default = NULL).
#' @param SE 	Standard error of estimate for each study.
#' @param Method A character string indicating which method is used to estimate the between-study variance \eqn{\tau^2} and its square root \eqn{\tau}.
#'
#' @return A named vector with elements:
#' \describe{
#' \item{R2M}{Estimated R2-based total mediation effect.}
#' \item{SE}{Asymptotic standard errors for R2M.}
#' \item{CI_width}{Half width of the 95% confidence interval.}
#' \item{CI_lower}{Lower bound of the 95% confidence interval.}
#' \item{CI_upper}{Upper bound of the 95% confidence interval.}
#' \item{p1}{Number of selected mediators in subsample 1.}
#' \item{p2}{Number of selected mediators in subsample 2.}
#' \item{R_YX}{Estimated R2 between X and Y.}
#' }
#'
#'
#'
#' @importFrom meta metagen
#'
#' @export
#' @examples
#' Study <- c(1,2,3,4,5)
#' Effects <- c(0.1,0.2,0.3,0.4,1.5)
#' SE <- c(0.01, 0.02, 0.01, 0.08, 0.1)
#' Method <- "Fixed"
#' MetaR2M::R2_Meta(Effects=Effects, Study=Study, SE=SE, Method=Method)
#' Method <- "DL"
#' MetaR2M::R2_Meta(Effects=Effects, Study=Study, SE=SE, Method=Method)
R2_Meta <- function(Effects, Study=NULL, SE, Method=c("Fixed", "REML", "PM", "DL", "ML")){
  if(is.null(Study)){
    study_names <- seq_len(length(Effects))
  }else{
    study_names <- Study
  }

  # Construct the data frame for the meta-analysis
  resultDF_S <- data.frame(Study = study_names, EstR2_S_tem = Effects, EstSD_S_tem = SE)


  if(Method == "Fixed"){
    Method_tem <- "DL"
  }else{
    Method_tem <- Method
  }

  maFixed <- meta::metagen(TE=resultDF_S$EstR2_S_tem, seTE=resultDF_S$EstSD_S_tem, studlab=resultDF_S$Study, method.tau = Method_tem)
  # maFixed <- mvmeta::mvmeta(resultDF_S$EstR2_S_tem ~ 1, S=resultDF_S$EstSD_S_tem^2, method=Method)

  if(Method == "Fixed"){
    meta_effect <- as.numeric(maFixed$TE.common)
    meta_CI_lower <- as.numeric(maFixed$lower.common)
    meta_CI_upper <- as.numeric(maFixed$upper.common)
    Q <- as.numeric(maFixed$Q)
    Q_pval <- as.numeric(maFixed$pval.Q)
    if(Q_pval <= 0.05){
      message("p value of Q is <= 0.05, try random-effects model")
    }


  }else{
    meta_effect <- as.numeric(maFixed$TE.random)
    meta_CI_lower <- as.numeric(maFixed$lower.random)
    meta_CI_upper <- as.numeric(maFixed$upper.random)
    Q <- as.numeric(maFixed$Q)
    Q_pval <- as.numeric(maFixed$pval.Q)
    if(Q_pval > 0.05){
      message("p value of Q is <> 0.05, try fixed-effects model")
    }
  }



  return(round(c(meta_effect=meta_effect, meta_CI_lower=meta_CI_lower,
                 meta_CI_upper=meta_CI_upper, Q=Q, Q_pval=Q_pval, Method=Method), 4))
}





