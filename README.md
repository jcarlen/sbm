# Code to fit various stochastic blockmodels 

## KLOptimization 
Contains standalone c++ code to fit degree-corrected (and not degree-corrected) stochastic blockmodels for a single (not multilayer) network.  KLOptimization.cpp is the original script from Karrer & Newman 2010 (http://www-personal.umich.edu/~mejn/dcsbm/) for undirected networks. KLOptimization_directed.cpp is my extension of that code to directed networks. See the comments of those files for implementation details.

## sbmt (TDD-SBM)

This **sbmt** folder contains an Rcpp package for stochastic block modeling to implement the TDD-SBM we introduced in ??. It builds off the KLOptimization code of Karrer and Newman, using a direct extension of the Kergnighan-Lin algorithm for multilyaer, specifically time-sliced, networks. However, the package is more general than fitting TDD-SBM, as it can also fit undirected networks, with or without degree correction, and with any type of non-negative edge weights. 

To install:
In terminal, from folder containing sbm > R CMD build sbm; R CMD check sbm; R CMD install sbm

This package is a work in progress and any suggestions for improvements or found bugs are appreciated. 
