# BiasedMatch 0.1.0

* Initial release implementing the two-urn biased sampling model of Puza and Bonfrer (2018).
* `DISTM()` returns the full distribution, mean and variance of the number of matches.
* `PROBM()` returns the probability of a given number of matches under biased without-replacement sampling.
* `Vmf()` gives the closed-form variance of the number of matches under unbiased sampling (`w = 1`).
* `inferW()` returns the maximum likelihood estimate, posterior mode and a p-value for the bias parameter `w` given an observed number of matches, reproducing Problems 1-3 of the paper.
* `likW()` returns standardised likelihood and posterior curves for the bias parameter.
* `MHALG()`, `RMMHALG()` and `SMRMMHALG()` provide Metropolis-Hastings samplers for hierarchical inference on the bias parameter.
* `restable()` summarises posterior draws and `estprop()` reports convergence and serial-correlation diagnostics.
