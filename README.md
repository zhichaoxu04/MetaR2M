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

- Download [MetaR2M](https://github.com/zhichaoxu04/MetaR2M) package from Github using:

<!-- -->

    git clone https://github.com/zhichaoxu04/crSKAT.git

- Or, install [MetaR2M](https://github.com/zhichaoxu04/MetaR2M) package in R directly

  - First, install [devtools](https://devtools.r-lib.org) in R from CRAN:
    ``` r
    install.packages("devtools")
    ```
  - Then, install [MetaR2M](https://github.com/zhichaoxu04/MetaR2M) using the `install_github` function and load the package:
    ``` r
    devtools::install_github("zhichaoxu04/crSKAT")
    library(crSKAT)
    ```
- Make sure that all the required packages have been installed or updated. Here are some of the required packages:
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

