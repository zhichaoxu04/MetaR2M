# MetaR2M: A novel meta-analysis framework for high-dimensional $R^2$-based mediation effect

# Introduction
## Meta-analysis
Meta-analysis is a statistical technique used to combine the results of multiple studies to arrive at a comprehensive understanding of a particular field or topic. This method is especially prevalent in fields like medicine, psychology, and social sciences, where individual studies might have varying outcomes or small sample sizes. By aggregating data from several studies, a meta-analysis can provide more robust conclusions, identify patterns, and offer insights that might not be apparent from individual studies. 
<div align="center"><img src="man/Figure/MetaAnalysis.png" ></div>
</br>

## Fixed/Random-effects model
Fixed-effects models require the assumption that the true effects of interest are identical across all studies or cohorts. Random-effects models are used when there is heterogeneity across the studies included in the meta-analysis. In meta-analysis, the choice between a fixed-effect and a random-effects model is fundamental and depends on the underlying assumptions about the nature of the effects being analyzed.

The fixed-effect model operates under the assumption that there is one true effect size that is common to all the studies being analyzed. In essence, it posits that any observed differences in effect sizes across studies are solely due to sampling error. This model is particularly appropriate when the meta-analysis includes studies that are highly similar in terms of participants, interventions, and outcomes, suggesting that they are all tapping into the same underlying effect. The key advantage of this model is its simplicity and increased statistical power due to the assumption of a single underlying effect. However, its major limitation is that it cannot account for variability beyond chance among different studies, which can lead to biased results if this assumption is violated.

On the other hand, the random-effects model acknowledges and accommodates heterogeneity in effect sizes across studies. This model assumes that the studies in the meta-analysis are estimating different, yet related, effect sizes. These differences could arise from variations in study populations, methodologies, or other contextual factors. The random-effects model includes both within-study sampling error and between-study variance in its calculations. This approach is more flexible and realistic in scenarios where study heterogeneity is expected. It provides a more generalized conclusion, applicable to a broader context beyond the specific studies included in the meta-analysis. The trade-off, however, is that this model often has less statistical power compared to the fixed-effect model due to the additional variance component that needs to be estimated.

<div align="center"><img src="man/Figure/FRmodels.png" ></div>
</br>

## $Q$ statistic for heterogeneity
The $Q$ statistic is a key measure in meta-analysis, primarily used to test for heterogeneity among the included studies. It is calculated as the weighted sum of squared differences between individual study effects and the overall effect estimate, where the weights are typically the inverse of each study's variance. The underlying principle of the $Q$ statistic is to determine whether the observed differences in study outcomes are greater than what would be expected by chance alone. Under the null hypothesis, which assumes homogeneity among the study effects (i.e., any differences are due to sampling error), the $Q$ statistic follows a chi-square distribution with degrees of freedom equal to the number of studies minus one. A significant $Q$ (usually indicated by a p-value less than a conventional threshold like 0.05) suggests the presence of heterogeneity, meaning that the variation in study results cannot be fully attributed to random chance. However,  $Q$'s power to detect heterogeneity depends on several factors, including the number of studies and their size. Due to this, $Q$ is often complemented with other measures, such as the $I^2$ statistic, which quantifies the proportion of total variation in study estimates due to heterogeneity.

# Get started
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
    devtools::install_github("zhichaoxu04/MetaR2M")
    library(MetaR2M)
    ```
- Make sure that all the required packages have been installed or updated. Here are some of the required packages:
  - [RsqMed](https://cran.r-project.org/web/packages/RsqMed/index.html): An implementation of calculating the R-squared measure as a total mediation effect size measure and its confidence interval for moderate- or high-dimensional mediator models. It gives an option to filter out non-mediators using variable selection methods. The original R package is directly related to the paper [Yang et al (2021)](https://pubmed.ncbi.nlm.nih.gov/34425752/). The new version contains a choice of using cross-fitting, which is computationally faster. The details of the cross-fitting method are available in our paper [Xu et al (2023)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9934518/).
  - [SIS](https://cran.r-project.org/web/packages/SIS/index.html): Variable selection techniques are essential tools for model selection and estimation in high-dimensional statistical models. Through this publicly available package, we provide a unified environment to carry out variable selection using iterative sure independence screening (SIS) ([Fan and Lv (2008)](https://academic.oup.com/jrsssb/article/70/5/849/7109492)) and all of its variants in generalized linear models ([Fan and Song (2009)](https://projecteuclid.org/journals/annals-of-statistics/volume-38/issue-6/Sure-independence-screening-in-generalized-linear-models-with-NP-dimensionality/10.1214/10-AOS798.full)) and the Cox proportional hazards model ([Fan, Feng and Wu (2010)](https://projecteuclid.org/ebooks/institute-of-mathematical-statistics-collections/Borrowing-Strength--Theory-Powering-Applications--A-Festschrift-for/chapter/High-dimensional-variable-selection-for-Coxs-proportional-hazards-model/10.1214/10-IMSCOLL606)).
  - [HDMT](https://cran.r-project.org/web/packages/HDMT/index.html): A multiple-testing procedure for high-dimensional mediation hypotheses. Mediation analysis is of rising interest in epidemiology and clinical trials. Among existing methods for mediation analyses, the popular joint significance (JS) test yields an overly conservative type I error rate and therefore low power. Methods used in the package refer to [James Y. Dai, Janet L. Stanford & Michael LeBlanc (2020)](https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1765785).
  - [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html): A Grammar of Data Manipulation: a fast, consistent tool for working with data frame like objects, both in memory and out of memory.
  - [meta](https://cran.r-project.org/web/packages/meta/index.html): User-friendly general package providing standard methods for meta-analysis and supporting Schwarzer, Carpenter, and Rücker, ["Meta-Analysis with R" (2015)](https://link.springer.com/book/10.1007/978-3-319-21416-0): - common effect and random effects meta-analysis; - several plots (forest, funnel, Galbraith / radial, L'Abbe, Baujat, bubble); - three-level meta-analysis model; - generalised linear mixed model; - Hartung-Knapp method for random effects model; - Kenward-Roger method for random effects model; - prediction interval; - statistical tests for funnel plot asymmetry; - trim-and-fill method to evaluate bias in meta-analysis; - meta-regression; - cumulative meta-analysis and leave-one-out meta-analysis; - import data from 'RevMan 5'; - produce forest plot summarising several (subgroup) meta-analyses.

## Toy Example

Let's suppose we aim to perform a meta-analysis using data from four distinct cohorts within a substantial genetic database. Our goal is to assess the mediating impact of gene expression on the age-related variations observed in systolic blood pressure (BP). For illustrative purposes, we will utilize the example data preloaded in our package as a representative sample for this analysis.

```r
library(MetaR2M)
```

