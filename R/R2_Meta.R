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
#' \item{meta_effect}{Estimated R2-based total mediation effect using meta-analysis.}
#' \item{meta_CI_lower}{Lower bound of the 95% confidence interval.}
#' \item{meta_CI_upper}{Upper bound of the 95% confidence interval.}
#' \item{Q}{Statistic used to test for heterogeneity among the studies included in the analysis.}
#' \item{Q_pval}{The p-value associated with the statistic.}
#' \item{Method}{Method is used to estimate the between-study variance in random effects model or fixed effects model.}
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
    meta_effect <- round(as.numeric(maFixed$TE.common), 4)
    meta_CI_lower <- round(as.numeric(maFixed$lower.common), 4)
    meta_CI_upper <- round(as.numeric(maFixed$upper.common), 4)
    Q <- round(as.numeric(maFixed$Q), 4)
    Q_pval <- round(as.numeric(maFixed$pval.Q), 4)
    if(Q_pval <= 0.05){
      message("p value of Q is <= 0.05, try random-effects model")
    }


  }else{
    meta_effect <- round(as.numeric(maFixed$TE.random), 4)
    meta_CI_lower <- round(as.numeric(maFixed$lower.random), 4)
    meta_CI_upper <- round(as.numeric(maFixed$upper.random), 4)
    Q <- round(as.numeric(maFixed$Q), 4)
    Q_pval <- round(as.numeric(maFixed$pval.Q), 4)
    if(Q_pval > 0.05){
      message("p value of Q is <> 0.05, try fixed-effects model")
    }
  }



  return(c(meta_effect=meta_effect, meta_CI_lower=meta_CI_lower,
                 meta_CI_upper=meta_CI_upper, Q=Q, Q_pval=Q_pval, Method=Method))
}





