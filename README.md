Binary AddiVortes
===========

An implementation of the Binary (Bayesian) Additive Voronoi Tessellation (Binary AddiVortes) algorithm in R.

Copyright (C) 2024

Adam Stone & John Paul Gosling

Department of Mathematical Sciences, Durham University
 
Setup Instructions
------------------

# Installing the Binary Algorithm

To install the AddiVortes functions in R, you can run the following code: 

```r

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)

source_url("https://raw.githubusercontent.com/Adam-Stone2/Binary_AddiVortes/main/Binary_AddiVortes.R")

```
The following function can then be used in Rstudio.

```r
 AddiVortes_Algorithm(y, x, m = 200, max_iter = 1200, burn_in = 200,
                       k = 3, sd = 0.8, Omega = 3, lambda_rate = 25,
                        YTest, XTest, IntialSigma = "Linear", thinning = 1)
```
### Arguments

`y`- Dependent variable for training (in sample) data. A numerical vector with index of output value corresponding to the row of observation in x.

`x`- Explanatory variables for training (in sample) data. A matrix with (as usual) rows corresponding to observations and columns to variables.

`m`- Number of tessellations.

`max_iter`- Number of iterations of the MCMC backfitting algorithm.

`Burn_in`- Number of iterations discarded before sampling posterior. 0 < burn_in < max_iter.

`k`- Number of prior standard deviations E(Y∣x)=f(x) is away from +/-.5. The response (y) is internally scaled to range from -.5 to .5. 

`sd` - Standard deviation of the normal distribution for new coordinates in tessellations.

`Omega`- Probability of covariate being included as a dimension in Tessellation prior. 

`lambda_rate`- The rate of the number of centres in a tessellation for Poisson distribution in tessellation prior.

`YTest`- Dependent variable for test (out of sample) data. These should have the same structure as y.

`XTest`- Explanatory variables for test (out of sample) data. These should have the same structure as x.

`IntialSigma`- Either "Linear" or "Naive". When “Naive”, a rough estimate of σ corresponds to the sample standard deviation of the transformed training response values.
If “Linear”, the rough estimate of σ is based on the residual standard deviation from a least-squares linear regression of Y on the original X variables.

`thinning`- Retaining every posterior post burn in sample equal to thinning. Default retains every posterior sample post burn in.

### Values

`In_sample_Accuracy`- The RMSE of the in-sample estimations.

`Out_of_sample_Accuracy`- The RMSE of the out-of-sample estimations.

Car Datasets
-----------------------------

To import the car dataset used in the paper and run Heteroscedastic AddiVortes in Rstudio one can run the following code:

```r
# Load necessary library if you haven't already
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
# Load the necessary library
library(readr)

Car_dataset <- read_csv("https://raw.githubusercontent.com/anonymous2738/Heteroscedastic_AddiVortes/main/Car_dataset.csv")
X_Car_dataset <- as.matrix(Car_dataset[,2:15])
Y_Car_dataset <- as.numeric(as.matrix(Car_dataset[,1]))


n <- length(Y_Car_dataset)
TrainSet <- sort(sample.int(n,4*n/5))
TestSet <- 1:n
TestSet <- TestSet[! TestSet %in% TrainSet]

AddiVortes_Algorithm(Y_Car_dataset[TrainSet],X_Car_dataset[TrainSet,],
                     200,40,2000,200,6,0.85,3,0.8,3,25,
                     Y_Car_dataset[TestSet],X_Car_dataset[TestSet,],
                     IntialSigma = "Linear")
```

Reproducing Figures in the paper 
---------------------------

**Warning:** These figure use parallel processing using up to 10 cores at a time, producing these figure is only recommended if you have 12+ cores.

To reproduce the figures 2-4 in the paper, source the following Github page by running:

```r

source_url("https://raw.githubusercontent.com/adam-stone2/Binary_AddiVortes/main/Figures2-4.R")

```
To reproduce the figures 5 in the paper, source the following Github page by running:

```r

source_url("https://raw.githubusercontent.com/anonymous2738/Heteroscedastic_AddiVortes/main/Figure5.R")

```
To reproduce the figures 6 in the paper, source the following Github page by running:

```r

source_url("https://raw.githubusercontent.com/anonymous2738/Heteroscedastic_AddiVortes/main/Figure6.R")

```
To reproduce the figures 7 in the paper, source the following Github page by running:

**Warning:** This Figure takes a long time to run (+3hours).

```r

source_url("https://raw.githubusercontent.com/anonymous2738/Heteroscedastic_AddiVortes/main/Figure7.R")

```
