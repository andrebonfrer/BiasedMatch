# BiasedMatch

`BiasedMatch` implements functions for biased matching problems where two
overlapping sets of items are sampled without replacement and items
common to both sets are sampled with higher probability.  The package
computes the probability mass function for the number of matching
items, its moments and provides Metropolis–Hastings algorithms to
perform Bayesian inference for the weight parameters.  The code is based
on the methods described by Puza and Bonfrer (2018).

## Installation

The package is distributed as a source package.  To install it from a
local directory, run the following in R:

```r
# install required dependencies first
install.packages(c("BiasedUrn", "mvtnorm", "nortest", "lawstat", "coda"))

# then install BiasedMatch from the local tar.gz
install.packages("path/to/BiasedMatch", repos = NULL, type = "source")
```

## Usage

The probability of obtaining a particular number of matches can be
computed with `PROBM()`, and the entire distribution with `DISTM()`.
Metropolis–Hastings algorithms for parameter inference are provided by
`MHALG()`, `RMMHALG()` and `SMRMMHALG()`.  Summaries and diagnostics
are available via `restable()` and `estprop()`.

```r
library(BiasedMatch)

# probability mass for 3 matches
PROBM(m = 3, N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 2)

# full distribution
DISTM(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 2)
```

See the package vignette for a detailed introduction and examples.

## Licence

This package is licensed under the MIT License (see `LICENSE`).