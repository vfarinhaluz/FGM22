# Risk classification numerical simulations

## Summary
This repository contains the code for numerical simulations (under active development) for the paper "Risk Classification in Insurance Markets with Risk and Preference Heterogeneity", 
by Vitor Farinha Luz, Piero Gottardi and Humberto Moreira.

Contact: [vitor.farinhaluz@ubc.ca](mailto:vitor.farinhaluz@ubc.ca)

Research assistance by Pedro Melgaré.

---
## Algorithm description:
Matlab price iteration:

---
## Files in repository:
BaselineEquilibrium.m: 
This code computes an approximate competitive equilibrium for the two-dimensional
model discussed in the main paper, allowing for any finite number of contracts and 
possible type realizations in each dimension. 

BaselineSignalAnalysis.m:
This code numerically simulates the effect of signal disclosure 
for the two-dimensional model discussed in the main paper, allowing for any 
finite number of contracts and possible type realizations in each dimension. 

SignalAnalysisGridRhoRisk.m:
This code creates a grid for two parameters of the model:
- delta: the variation in risk aversion
- deltaRisk: the variation in risk levels

For each pair of deltaRisk and delta, the code simulates the effect of signal
disclosure and stores all the simulated data into a large Signal Simulation
array object.

This code also generates some heatmaps to display variables of interest,
such as 1.  1. the error term in the iteration algorithm, 2. the mass of 
consumers with interim improvement, and 3. the aggregate welfare effect of 
the signal disclousure

SignalAnalysisGridRhoDelta.m:
Same than SignalAnalysisGridRhoRisk.m, but this code considers a grid of rho0 (the mean value for risk aversion) and delta.


\src\:

HyperParameters.m: 
defines HyperParameters class, which collects parameters for the iteration algorithm

pdfEinav.m:
implements the pdf function for a lognormal distribution with mean and covariance parameters set to Einav et al (2013) paper's values.

pdfEinavFlex.m:
does the same than pdfEinav, but allows the user to input the first moment of the distribution.

Model.m:
1. defines Model class, which collects parameters for the model to be solved
2. implements several functions to compute consumer demand, average costs, consumer's utility levels, and others

Solver.m:
defines Solver class, which collects a Model and HyperParameters. Also in this file is the iteration algorithm used for finding equilibria. 

SimulationResult.m:
1. defines SimulationResult class, which collects equilibrium objects (prices, utilities) and the error (ie, the maximum |price-AC| over all contracts) across iterations,
as well as model parameters.
2. implements functions to plot equilibrium objects, such as prices, demand, and type assignment

SignalSimulation.m:
This file defines the class SignalSimulation, solves a binary signal disclousure exercise, compute welfare effects, plots changes in equilibrium objects such as prices and utility levels.
