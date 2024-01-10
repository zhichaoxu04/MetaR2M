#' @title Simulated data for tutorial
#'
#' @description simulated 4 independent studies for meta-analysis
#'
#' @format A list of 4 studies with different sample sizes and each of them contains:
#' \describe{
#'   \item{X}{The exposure mimics age}
#'   \item{Y}{The outcome mimics systolic blood pressure}
#'   \item{M}{The mediators}
#' }
#' @references <https://www.biorxiv.org/content/10.1101/2023.02.06.527391v1.abstract>
"exampleData"
#
# N <- 2000
# Outcome <- rnorm(N, 100, sd = 15)
# Exposure <- rnorm(N, 60, sd = 10)
# M <- matrix(0, nrow = N, ncol = 50)
# for (i in 1:ncol(M)){M[, i] <- rnorm(N, 0, 1) * Exposure + rnorm(N, 0, 1)}
# colnames(M) <- paste0("M_", 1:50)
# exampleData <- list(study1 = cbind(data.frame(Y = Outcome[1: 300], X = Exposure[1: 300]), M[1:300, ]),
#                     study2 = cbind(data.frame(Y = Outcome[301: 800], X = Exposure[301: 800]), M[301: 800, ]),
#                     study3 = cbind(data.frame(Y = Outcome[801: 1200], X = Exposure[801: 1200]), M[801: 1200, ]),
#                     study4 = cbind(data.frame(Y = Outcome[1201: 2000], X = Exposure[1201: 2000]), M[1201: 2000, ]))
# usethis::use_data(exampleData, overwrite = TRUE)
