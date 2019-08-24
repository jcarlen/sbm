# Code to fit time-dependent stochastic block models 

## KLOptimization 
Contains standalone c++ code to fit degree-corrected (and not degree-corrected) stochastic block models for a single-layer (not time-sliced) network. We include it here because our code to fit time-dependent discrete-membership SBM (in sbmt) build off this code. In KLOptimization, KLOptimization.cpp is the original script from Karrer & Newman 2010 (http://www-personal.umich.edu/~mejn/dcsbm/) for undirected networks. KLOptimization_directed.cpp is my extension of that code to directed networks. See the comments of those files for implementation details.

## sbmt (TDD-SBM)

This **sbmt** folder contains an Rcpp package to implement the time-dependent discrete stochastic block model (TDD-SBM) we introduce in ??. It builds off the KLOptimization code of Karrer and Newman, using a direct extension of the Kergnighan-Lin algorithm for multilayer (specifically time-sliced) networks. However, the package is more general than fitting TDD-SBM, as it can also fit undirected networks, with or without degree correction, with or without multilpe layers (time slice), and with any type of non-negative edge weights.

To install:
In terminal, from folder containing sbmt > R CMD build sbmt; R CMD check sbmt; R CMD install sbmt

This package is a work in progress and any suggestions for improvements or found bugs are appreciated. 
