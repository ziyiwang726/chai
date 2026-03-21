# chai
Conditional Hypothesis testing using Auxiliary Information

[![R-CMD-check](https://img.shields.io/badge/R-CMD--check-passing.svg)](https://github.com/your_username/chai/actions)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
## Overview

`chai` is an R package developed for a covariate-informed statistical framework. It leverages auxiliary information to enhance the statistical power of multiple hypothesis testing while controlling the false discovery rate (FDR) of high-dimensional data (such as 16S rRNA and WGS microbiome sequencing).

## Authors

**Ziyi Wang, Satabdi Saha, Christine B. Peterson, Yushu Shi**


## Installation

You can install the development version of `chai` from GitHub using `remotes`:

```R
# Install remotes if you haven't already
if (!require("remotes")) install.packages("remotes")

# Install chai
remotes::install_github("ziyiwang726/chai")
```

## Example

### Generate a simulation
We assumed a total 1000 hypotheses, where the first 950 are null hypotheses with $z \sim N(0, 1)$, and the remaining 50 are alternatives with $z \sim N(3, 1)$. 
Correspondingly, $x$ was used as side information: for the first 950 hypotheses, $x \sim N(3, 1)$, and for the last 50 hypotheses, $x \sim N(6, 1).
Thus, the true alternative hypothesis indices are from 951 to 1000. 

```R
set.seed(123)
n = 1000; n0 <- 950; n1 <- 50
z0 <- rnorm(n0, mean = 0, sd = 1)
x0 <- rnorm(n0, mean = 3, sd = 1)
z1 <- rnorm(n1, mean = 3, sd = 1)
x1 <- rnorm(n1, mean = 6, sd = 1)
z <- c(z0, z1)
X <- c(x0, x1)
gt <- seq((n0+1), n)    # Ground Truth
```

We then fit both $z$ and $x$ into the funciton `chai`:

```R
# Fit the model
res <- chai(z, X, K_vec = 2:6, B = 100)
```

Examined which hypotheses were rejected at $q = 0.05$ (target FDR level) using the function `clfdrselect`:
```R
# Check the rejections
clfdrselect(res$clFDR, q = 0.05)
```

Based on these rejections and the true ground truth, we evaluated the model's FDP and statistical power:
```R
# Check performance with ground truth
performance(gt, clfdrselect(res$clFDR, q = 0.05))
```

