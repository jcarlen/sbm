# This script illustrates how sbmt's KL algorithm can be easily parallelized
# NOTE: This works on my mac but hasn't been tested on other OS so is not implemented in the package
# ---------------------------------------------------------------------------------------------------
# libraries ----

#devtools::install_github("jcarlen/sbm", subdir = "sbmt")

library(sbmt)
library(dplyr)
library(parallel)

# load sample data from package ----

data("la_byhour_edgelist")

# params ----

kl = 10 #klPerNetwork
if (kl < 1) {kl = 1; cat("invalid kl, using 1 kl run")}
ncluster = 2
if (ncluster > detectCores()) {ncluster = 1; cat("invalid ncluster, using 1")}
if (!exists(".Random.seed")) set.seed(NULL)
seed = 1

# setup cluster ----

cl <- makeCluster(ncluster, type='FORK')
clusterSetRNGStream(cl, seed)

# setup cluster runs (each cluster does kl/ncluster runs) ----

sbmt.chunk <- function(x) {sbmt(edgelist.time = la_byhour_edgelist, maxComms = 2, degreeCorrect = 3, directed = T, 
                                klPerNetwork = ceiling(kl/ncluster), tolerance = 1e-4, seed = NULL, seedComms = NULL)}

# for sbm function (static SBM) the chunked input would just look like:
# sbm.chunk <- function(x) {sbm(edgelist = la_byhour_edgelist[[18]], maxComms = 2, degreeCorrect = 0, directed = T, klPerNetwork = ceiling(kl/ncluster), tolerance = 1e-4, seed = NULL, seedComms = NULL)}

# Do all KL runs and select best run for final result ----

sbmt.chunk.out = parLapply(cl, seq_len(ncluster), sbmt.chunk)

stopCluster(cl)

best.run = which.max(sapply(sbmt.chunk.out, "[[", "llik"))
Results = sbmt.chunk.out[[best.run]]