# Risk classification numerical simulations

## Summary
This repository contains the code for numerical simulations (under active development) for the paper "Risk Classification in Insurance Markets with Risk and Preference Heterogeneity", 
by Vitor Farinha Luz, Piero Gottardi and Humberto Moreira.

Contact: [vitor.farinhaluz@ubc.ca](mailto:vitor.farinhaluz@ubc.ca)

Matlab research assistance by Pedro Melgaré.

This code can also be accessed at:
[https://github.com/vfarinhaluz/FGM22](https://github.com/vfarinhaluz/FGM22)

---

## Algorithm description:

An equilibrium in the model is constituted of a price function **p()** and an allocation which maps from the set of types to coverage levels (the paper often uses the inverse mapping, which is a type assignment function).
A type contains both risk (mu) and risk aversion (rho). The algorithm uses a finite grid as the set of possible types.

The algorithm starts from a initial guessed price function and iterates a properly defined operator until convergence is achieved, which represents that an approximate equilibrium has been found.

More precisely, for a given price function p_t, the algorithm calculates the optimal coverage level for each type and defines a new price function  p' as follows:

1. For each contract traded under p_t, the new price p' is updated towards  the average risk of agents purchasing it. 
2. The price for each non-traded contract x' is also updated,  but on the basis of the riskiness of the type most willing to pay for it.

The new price p_(t+1) is given by 

p_(t+1) = p_t + alpha (p' - p_t),

for some parameter alpha in (0,1) (this parameter is chosen to make sure prices don't fluctuate too much in each step but still converge).

The simulation results are reported in Section 5.5 of the main paper, and Section 1 of the online appendix.

---
## Files and folders in repository

### BaselineEquilibrium.m: 
This code computes an approximate competitive equilibrium for the two-dimensional
model discussed in the main paper, allowing for any finite number of contracts and 
possible type realizations in each dimension. 

### BaselineSignalAnalysis.m:
This code numerically simulates the effect of signal disclosure 
for the two-dimensional model discussed in the main paper, allowing for any 
finite number of contracts and possible type realizations in each dimension. 


### SignalAnalysisGridRhoDelta.m:

This code creates a grid varying the following parameters of the model:
- delta: the variation in risk aversion
- rho_0 : the mid-point of the risk aversion range

For each pair of delta and rho_0, the code simulates the effect of signal
disclosure and stores all the simulated data into a large Signal Simulation
array object.

This code also generates heatmaps to display variables of interest,
such as: the error term in the iteration algorithm, the mass of 
consumers with interim improvement, and the aggregate welfare effect of 
the signal disclosure.

### SignalAnalysisGridRhoRisk.m:

This code creates a grid varying the following:
- delta: the variation in risk aversion
- (mu_H - mu_L) : the variation in risk levels

For each pair of deltaRisk and delta, the code simulates the effect of signal
disclosure and stores all the simulated data into a large Signal Simulation
array object.

This code also generates heatmaps to display variables of interest.

### \src\:

Folder containing the following:

1. HyperParameters.m: 
defines HyperParameters class, which collects parameters for the iteration algorithm

1.  pdfEinavFlex.m: implements the pdf function for a lognormal distribution with covariance parameters set to Einav et al (2013) paper's values, allowing the user to input the first moment of the distribution.

2.  Model.m: 
    1.  defines Model class, which collects parameters for the model to be solved
    2. implements several functions to compute consumer demand, average costs, consumer's utility levels, and others

3.  Solver.m:
defines Solver class, which collects a Model and HyperParameters. Also in this file is the iteration algorithm used for finding equilibria. 

1.  SimulationResult.m:
    1. defines SimulationResult class, which collects equilibrium objects (prices, utilities) and the error (ie, the maximum |price-AC| over all contracts) across iterations, as well as model parameters.
    2. implements functions to plot equilibrium objects, such as prices, demand, and type assignment

2.  SignalSimulation.m:
This file defines the class SignalSimulation, solves a binary signal disclosure exercise, compute welfare effects, plots changes in equilibrium objects such as prices and utility levels.
