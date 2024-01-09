#' R2_Meta.R
#'
#' A novel meta-analysis framework of the R-squared (R2)-based mediation effect estimation for high-dimensional omics mediators.
#'
#' @param Effects Estimated R2M for each study.
#' @param Study Study name (Default = NULL).
#' @param SE Asymptotic standard errors for each study.
#' @param Method Methods for the meta-analysis.
#'
#' @return A named vector with elements:
#'     \describe{
#'       \item{R2M}{Estimated R2-based total mediation effect.}
#'       \item{SE}{Asymptotic standard errors for R2M.}
#'       \item{CI_width}{Half width of the 95% confidence interval.}
#'       \item{CI_lower}{Lower bound of the 95% confidence interval.}
#'       \item{CI_upper}{Upper bound of the 95% confidence interval.}
#'       \item{p1}{Number of selected mediators in subsample 1.}
#'       \item{p2}{Number of selected mediators in subsample 2.}
#'       \item{R_YX}{Estimated R2 between X and Y.}
#'       \item{R_YM}{Estimated R2 between M and Y.}
#'       \item{R_YXM}{Estimated R2 between XM and Y.}
#'       \item{SOS}{Shared Over Simple measure, which is R2M/R_YX.}
#'       \item{SOS_CI_lower}{Lower bound of SOS.}
#'       \item{SOS_CI_upper}{Upper bound of SOS.}
#'       \item{Time}{Computational time (in minutes).}
#'       \item{SampleSize}{Sample size used in the analysis}
#'       \item{NumMeds}{Number of total mediators for the input.}
#'     }
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
#' Method <- "fixed"
#' MetaR2M::R2_Meta(Effects=Effects, Study=Study, SE=SE, Method=Method)
#' Method <- "REML"
#' MetaR2M::R2_Meta(Effects=Effects, Study=Study, SE=SE, Method=Method)
R2_Meta <- function(Effects, Study=NULL, SE, Method=c("Fixed", "REML", "PM", "DL", "ML")){
  if(is.null(Study)){
    study_names <- seq_len(length(Effects))
  }else{
    study_names <- Study
  }

  # Construct the data frame for the meta-analysis
  resultDF_S <- data.frame(Study = study_names, EstR2_S_tem = Effects, EstSD_S_tem = SE)
  maFixed <- meta::metagen(TE=resultDF_S$EstR2_S_tem, seTE=resultDF_S$EstSD_S_tem, studlab=resultDF_S$Study, method.tau = Method)
  # maFixed <- mvmeta::mvmeta(resultDF_S$EstR2_S_tem ~ 1, S=resultDF_S$EstSD_S_tem^2, method=Method)

  # if()
  maFixed$TE.common

  meta_effect <- as.numeric(maFixed$coefficients)

  return(c(meta_effect=meta_effect))
}





