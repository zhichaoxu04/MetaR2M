# MetaR2M: A novel meta-analysis framework for high-dimensional $R^2$-based mediation effect

## Introduction
### Meta-analysis
Meta-analysis is a statistical technique used to combine the results of multiple studies to arrive at a comprehensive understanding of a particular field or topic. This method is especially prevalent in fields like medicine, psychology, and social sciences, where individual studies might have varying outcomes or small sample sizes. By aggregating data from several studies, a meta-analysis can provide more robust conclusions, identify patterns, and offer insights that might not be apparent from individual studies. 
<div align="center"><img src="man/Figure/MetaAnalysis.png" ></div>
</br>

### Fixed/Random-effects model
Fixed-effects models require the assumption that the true effects of interest are identical across all studies or cohorts. Random-effects models are used when there is heterogeneity across the studies included in the meta-analysis.

<div align="center"><img src="man/Figure/FRmodels.png" ></div>
</br>

## Get started
Download and install following required R packages:

- Download [crSKAT](https://github.com/zhichaoxu04/crSKAT) package from
  Github using:

<!-- -->

    git clone https://github.com/zhichaoxu04/crSKAT.git

- Or, install [crSKAT](https://github.com/zhichaoxu04/crSKAT) package in
  R directly

  - First, install [devtools](https://devtools.r-lib.org) in R from
    CRAN:

    ``` r
    install.packages("devtools")
    ```

  - Then, install [crSKAT](https://github.com/zhichaoxu04/crSKAT) using
    the `install_github` function and load the package:

    ``` r
    devtools::install_github("zhichaoxu04/crSKAT")
    library(crSKAT)
    ```

- Make sure that all the required packages have been installed or
  updated. Here are some of the required packages:

  - [CompQuadForm](https://cran.r-project.org/web/packages/CompQuadForm/index.html):
    Compute the distribution function of quadratic forms in normal
    variables using Imhof’s method, Davies’s algorithm, Farebrother’s
    algorithm or Liu et al.’s algorithm.
  - [nleqslv](https://cran.r-project.org/web/packages/nleqslv/index.html):
    Solves a system of nonlinear equations using a Broyden or a Newton
    method with a choice of global strategies such as line search and
    trust region.
  - [ICSKAT](https://cran.r-project.org/web/packages/ICSKAT/index.html):
    Implements the Interval-Censored Sequence Kernel Association
    (ICSKAT) test for testing the association between interval-censored
    time-to-event outcomes and groups of single nucleotide polymorphisms
    (SNPs).
  - [bindata](https://cran.r-project.org/web/packages/bindata/index.html):
    Generates correlated artificial binary data.

## Toy Example

Let’s assume our objective is to study the possible association between
a specific gene and the time taken for a fracture to occur. In this
context, death acts as the competing risk. To delve into this premise,
we’ll simulate event times for 5,000 subjects based on a proportional
hazards model. Observations will be scheduled at four distinct times: 4,
11, 18, and 25. Each participant’s precise observation time will be
determined by a random draw from a Uniform(-0.25, 0.25) distribution
centered around these predefined times. Furthermore, there’s a 10%
chance for any participant to miss a scheduled visit.

The genetic dataset will include information on 50 single nucleotide
polymorphisms (SNPs) related to the gene in question. For every patient,
their minor allele count, which can be 0, 1, or 2, is documented for
each of these 50 SNPs. Alongside the genetic data, we also possess
details on non-genetic covariates: one being categorical and the other
one continuous.

``` r
library(crSKAT)
# ---- Initialization of random seed, sample size, the number of SNPs, and parameters
set.seed(1)
n <- 5000
q <- 50
alpha1 <- -0.058
alpha2 <- -0.035
beta1 <- 0.03
beta2 <- log(1 - exp(beta1 / alpha1)) * alpha2

# ---- Generate dataset
# Assume all SNPs have MAF of 0.2 in this toy example
gMat <- matrix(data=rbinom(n=n*q, size=2, prob=0.2), nrow=n)
gSummed <- matrix(data=apply(gMat, 1, sum), ncol=1)
# Two covariates (one categorical and one continuous)
xMat <- matrix(data=c(rbinom(n=n, size=1, prob=0.5), 
                      runif(n=n, min=0, max=10)), 
               nrow=n)
# Observation time windows
obsTimes <- seq(4, 28, 7)
# Generate data
outcomeDat <- crSKAT::genData(seed=2023, n=n, 
                              alpha1=alpha1, alpha2=alpha2, 
                              beta1=beta1, beta2=beta2,
                              obsTimes=obsTimes, probMiss=0.1, 
                              windowHalf=.25)
# Left/right exact time for each subject
lt <- outcomeDat$leftTimes
rt <- outcomeDat$rightTimes
# Indicator of right-censored or not
obsInd <- outcomeDat$deltaVecSimple
# Indicator of competing outcome or primary outcome
deltaVec <- outcomeDat$deltaVec

# ---- Perform inference
# make design matrix with cubic spline terms using one knot
dmats <- crSKAT::makeICdmat(xMat=xMat, lt=lt, rt=rt, 
                            obsInd=obsInd, quant_r=NULL, nKnots=1)
# Fit null model using the genotype information 
# only need to do this once for each genetic set 
# note: there is only gSummed on the SNPs used here, which will be constant
nullFit <- crSKAT::crSKAT_fit_null(init_beta=c(rep(0,10),.001), 
                                   leftDmat=dmats$left_dmat, rightDmat=dmats$right_dmat,
                                   deltaVec=deltaVec, leftTimes=lt, gSummed=gSummed, 
                                   allowSingular=TRUE, method="Broyden")

# Perform the crSKAT and crBurden test
crICSKATOut <- crSKAT::crSKAT(leftDmat=dmats$left_dmat, rightDmat=dmats$right_dmat, 
                              leftTimes=lt,deltaVec=deltaVec, gMat=gMat, 
                              gSummed=gSummed, null_beta=nullFit$beta_fit, 
                              pvalue=TRUE)
# p-value of crSKAT
crICSKATOut$p_SKAT
#> [1] 0.6765122
# p-value of crBurden test
crICSKATOut$p_burden
#> [1] 0.5926542

# --- Test another genotype matrix, we DO NOT need to fit the null model again
# Generate another set of SNPs
gMat <- matrix(data=rbinom(n=n*q, size=2, prob=0.1), nrow=n)
gSummed <- matrix(data=apply(gMat, 1, sum), ncol=1)
# Perform the crSKAT and crBurden test again
crICSKATOut <- crSKAT::crSKAT(leftDmat=dmats$left_dmat, rightDmat=dmats$right_dmat, 
                              leftTimes=lt,deltaVec=deltaVec, gMat=gMat, 
                              gSummed=gSummed, null_beta=nullFit$beta_fit, 
                              pvalue=TRUE)
# p-value of crSKAT
crICSKATOut$p_SKAT
#> [1] 0.9850078
# p-value of crBurden test
crICSKATOut$p_burden
#> [1] 0.277369
```
