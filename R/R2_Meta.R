#' R2_Meta.R
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
#' @return A list with the three elements:
#' \itemize{
#'   \item A named vector with elements:
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
#'   \item \code{m1}{Selected mediators in subsample 1.}
#'   \item \code{m2}{Selected mediators in subsample 2.}
#' }
#'
#'
#'
#' @importFrom stats rnorm residuals lm var qnorm cov
#' @importFrom SIS SIS
#' @importFrom dplyr select
#' @importFrom HDMT null_estimation fdr_est
#'
#' @export
#' @examples
#' p0 <- 20
R2_Meta <- function(Y, M, Covar, X, iter.max=3, nsis=NULL, seed=2024, FDR=FALSE, FDRCutoff=0.2, method="iSIS"){
  x <- 1
  return(x)
}





