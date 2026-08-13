# run_promo_MHalg.R
# Test promotions for MHalg, focusing only on the most interesting category: the product categories themselves
# To be run from this working directory.

rm(list=ls())
source("MHalg.R")               			# this will load up the scripts needed, including supporting libraries and scripts.

load("Ddt.RData")                			# saved data from earlier

# Length of MCMC chain
warmup <- 1000
numiter <- 1000

# Set up vectors needed for algorithm
# 1 = Retailer 1, 2 = Retailer 2
N1vec <- Ddt$N1								# assortment 1
N2vec <- Ddt$N2								# assortment 2
N0vec <- Ddt$N								# overlap of 1 and 2
M1vec <- Ddt$m1								# set chosen from assortment 1
M2vec <- Ddt$m2								# set chosen from assortment 2 
mvec  <- Ddt$m								# outcome: number of matches 

# set up categories 
categories <- unique(Ddt$category)
Cvec <- match(Ddt$category,categories)

t1 <- 0.7
t2 <- 0.5
t3 <- 0.3
t4 <- .7
t5 <- .7
t6 <- .7
t7 <- 0.6
t8 <- .9
t9 <- 0.9
t10 <- .9

tvector <- c(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)

set.seed(561); date()	#

source('initvector.R')
res=MHALG(B=warmup, J=numiter,  N1vec=N1vec,  N2vec=N2vec, N0vec=N0vec, M1vec=M1vec, M2vec=M2vec, mvec=mvec, Cvec=Cvec, eta=0.0001, tau=0.0001, mu0=0, sig0=10000, muinit=0, rinit=initvector,
del=tvector)

initvector <-  tail(res$rm,1)            # last state of chain

dump(c('res','initvector'), file= "initvector.R")
date()	# Took this long
res$rar  # acceptance rate

# plot results for visual inspection
matplot(res$rm[-(1:warmup),],type='l')
