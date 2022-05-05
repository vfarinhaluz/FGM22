%{
Project: Risk Classification in Insurance Markets with Risk and Preference
Heterogeneity, by Vitor Farinha Luz, Humberto Moreira, Piero Gottardi
Matlab research assistance: Pedro Melgare
April 2nd, 2022

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
%}

%% Setup
% Code to go over grid of parameters and check signal effect.
addpath('./src');
clear;close all; clc;

%% Baseline parameters:
nContracts = 250;
nRisk = 500;
nRiskAv = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Type distribution

% Average type levels:
rho_0= 8;
mu_0 = 5;

% Type distribution - Einav et al (2013) - Lognormal
distributionName='Einav et al (2013)';
distribution_pdf = @(mu, rho) pdfEinavFlex(mu, rho, mu_0, rho_0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation parameters
BehavMass = 0;
tolerance = 1e-4;
maxIterations = 1e4;
alpha = 0.001; % alpha defines the size of the step
hyperParameters = HyperParameters;
hyperParameters.Tolerance = tolerance;
hyperParameters.MaxIterations = maxIterations;
hyperParameters.Alpha = alpha;

%Creating grids 
grid_Delta_percent=linspace(.2,.9,20);
grid_DeltaRiskPercent= linspace(.2,.9,20);

%Create array varying all options. 
% Array dimensions are:
% lines: rho_0 values
% second (column): delta values
[Array_Delta_percent, Array_DeltaRiskPercent] = meshgrid(grid_Delta_percent, grid_DeltaRiskPercent);

%% Running grid
signalSimuArray=SignalSimulation.empty; %Empty array to store all simulations
length_grid=numel(Array_Delta_percent);

for i=1:length_grid
    % Risk aversion type support:
    delta= Array_Delta_percent(i) * rho_0;
    rhoL = rho_0-delta/2;
    rhoH = rho_0+delta/2;

    % Risk type support:
    deltaRisk = Array_DeltaRiskPercent(i) * mu_0;
    muH = mu_0 + deltaRisk/2 ;
    muL = mu_0 - deltaRisk/2 ;

    % MONOTONIC SIGNAL:
    signal_A_Fun=@(m) (m-muL)./(muH-muL); %Must be a vectorized function in [0,1]
    
%     % NON MONOTONIC SIGNAL
%     signal_A_Fun=@(m) SignalSimulation.nonMonotSignal(muL,muH,m,.3,.7); %Must be a vectorized function in [0,1]

    model=Model(muL, muH, rhoL, rhoH, nRisk, nRiskAv, nContracts, distribution_pdf, distributionName, BehavMass);
    solver = Solver(model,hyperParameters);
    signalSimuArray(i)= SignalSimulation.simulateSignal(signal_A_Fun,solver,"min");
    sprintf('#######  PERCENT COMPLETE: %.5g  ######',100*i/length_grid)
end

SignalSimulation.saveNoHistory(signalSimuArray, 'SignalSimulationGridDeltaRisk');

%% Creating a table object with all relevant information:
 
plotHeatMapsRiskDelta(signalSimuArray,1);
