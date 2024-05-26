# lmoments_redux

This repository stores replication material for the paper "Inference in parametric models with many L-moments". It is structured as follows:

* __methods__: contains generic functions to implement the main methods in the paper. These function are used by the replication material in the other folders.

* __montecarlo__: materials for replicating the Monte Carlo exercises in the paper. The scripts _gev.R_ and _gpd.R_ replicate, respectively, exercises using a Generalised Exponential and Generalised Pareto distribution. The subfolder __aux__ contains auxiliary scripts to these exercises that are used by the _gev.R_ and _gpd.R_ routines.

* __application__: materials for replicating the empirical application in the paper. The chunk _run_application.R_ replicates the results in the paper. The file _dataset_ridesharing.csv_ stores the data used in the application. The file _specifications_application.R_ is an auxiliary script that is called by _run_application.R_ and contains functions representing the specifications adopted in the application.