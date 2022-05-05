# Risk classification project
This repository contains the code for numerical simulations (under active development) for the paper "Risk Classification in Insurance Markets with Risk and Preference Heterogeneity", 
by Vitor Farinha Luz, Piero Gottardi and Humberto Moreira.

Code also being developed by:
- Pedro Melgaré: [pedro_melgare@hotmail.com](mailto:pedro_melgare@hotmail.com)
- Goh Yeow Chong: [gohyc1993@gmail.com](mailto:gohyc1993@gmail.com)

List of changes/outstanding issues
1. Removing OUTDATEDAVerageCost from Model.m
2. I do not believe demandArray is being used at all at the moment; should we delete it?
3. Removing UtilityMatrix, since CARA_UtilityMatrix computes the same object
4. plotCompareAverageCost were made for debugging purposes, I don't think we should include it in the final version
5. simulateSignal moved from Solver.m to SignalSimulation.m
6. updateModelSignal moved from Model.m to SignalSimulation.m
7. removed QuasiLin_UtilCompareSignal from SignalSimulation.m since it is not being used anywhere
8. removed demandSize and errorPricePercent from Model.m since they are not being used anywhere
9. in name of avoiding redundancies, I merged all plotting functions in SignalSimulation.m

Matlab price iteration:

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
