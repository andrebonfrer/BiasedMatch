# BiasedMatch

<!-- badges: start -->
<!-- badges: end -->

`BiasedMatch` implements the two-urn biased sampling model of Puza and Bonfrer
(2018). Two overlapping sets of items are sampled without replacement, with
items common to both sets given a selection weight that favours their
inclusion. The package computes the distribution of the number of matching
items, performs maximum likelihood and Bayesian inference for the bias
parameter, and provides Metropolis-Hastings samplers for hierarchical
extensions.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("andrebonfrer/BiasedMatch")
```

## Usage

The full distribution, mean and variance of the number of matches:

```r
library(BiasedMatch)

# Probability of exactly 3 matches (Example 1 of the paper, w = 2)
PROBM(m = 3, N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 2)

# Full distribution with mean and variance
DISTM(N1 = 37, N2 = 45, N = 16, m1 = 12, m2 = 8, w = 2)
```

Inference on the selection bias `w` from an observed number of matches
(reproducing the marketing application, Problem 3):

```r
inferW(m = 13, N1 = 292, N2 = 350, N = 192, m1 = 90, m2 = 59)
#> $mle              1.5265
#> $posterior_mode   1.3029
#> $p_value          0.1737
```

`likW()` returns the standardised likelihood and posterior curves for
plotting, and the `MHALG()` family of functions performs hierarchical
Metropolis-Hastings inference when several experiments are pooled.

See `vignette("biasedmatch")` for a full walkthrough that reproduces the
tables and figures of the paper.

## Reference

Puza, B. and Bonfrer, A. (2018). A series of two-urn biased sampling problems.
*Communications in Statistics - Theory and Methods*, 47(1), 80-91.
<https://doi.org/10.1080/03610926.2017.1300282>

## Licence

MIT (see the `LICENSE` file).
