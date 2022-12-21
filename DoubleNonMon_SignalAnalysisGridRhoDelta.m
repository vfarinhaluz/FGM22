%{
Project: Risk Classification in Insurance Markets with Risk and Preference
Heterogeneity, by Vitor Farinha Luz, Humberto Moreira, Piero Gottardi
Matlab research assistance: Pedro Melgare
April 2nd, 2022

This code creates a grid for two parameters of the model:
- rho_0: the midpoint of risk aversion
- delta: the variation in risk aversion

For each pair of rho_0 and delta, the code simulates the effect of signal
disclosure and stores all the simulated data into a large Signal Simulation
array object.
%}

%% Setup
% Code to go over grid of parameters and check signal effect.
addpath('./src');
clear;close all; clc;

%% Baseline parameters:
nContracts = 250;
nRisk = 30;
nRiskAv = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Type distribution

% Risk type support:
mu_0 = 5;
delta_mu = 0.8;
muH = mu_0 * (1+delta_mu/2);
muL = mu_0 * (1-delta_mu/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BehavMass = 0;
tolerance = 1e-4;
maxIterations = 1e4;
alpha = 0.001; % alpha defines the size of the step
hyperParameters = HyperParameters;
hyperParameters.Tolerance = tolerance;
hyperParameters.MaxIterations = maxIterations;
hyperParameters.Alpha = alpha;

% % MONOTONIC SIGNAL:
% signal_A_Fun=@(m) (m-muL)./(muH-muL); %Must be a vectorized function in [0,1]

% NON MONOTONIC SIGNAL
signal_A_Fun=@(m) SignalSimulation.nonMonotSignal(muL,muH,m,.3,.7); %Must be a vectorized function in [0,1]

%Creating grids
grid_Delta_percent=linspace(.2,.9,20);
grid_rho0= linspace(2,15,20);

%Create array varying all options. 
% Array dimensions are:
% lines: rho_0 values
% second (column): delta values
[Array_Delta_percent, Array_rho0] = meshgrid(grid_Delta_percent, grid_rho0);

%% Generating simulations over the parameter grid
length_grid=numel(Array_rho0);
signalSimuArray = SignalSimulation.empty;
for i=1:length_grid
    rho_0=Array_rho0(i);
    delta=Array_Delta_percent(i).*rho_0;
    rhoL = rho_0-delta/2;
    rhoH = rho_0+delta/2;

    % Type distribution - Einav et al (2013) - Lognormal
    distributionName='Einav et al (2013)';
    distribution_pdf = @(mu, rho) pdfEinavFlex(mu, rho, mu_0, rho_0);

    model=Model(muL, muH, rhoL, rhoH, nRisk, nRiskAv, nContracts, distribution_pdf, distributionName, BehavMass);
    solver = Solver(model,hyperParameters);
    signalSimuArray(i)= SignalSimulation.simulateSignal(signal_A_Fun,solver,"min");
    sprintf('#######  PERCENT COMPLETE: %.5g  ######',100*i/length_grid)
end

% Activate this line to save data on simulations:
% SignalSimulation.saveNoHistory(signalSimuArray, 'SignalSimulationGridRho0Delta');

%% Creates heatmaps to illustrate how parameters affect signal effect:
 
plotHeatMapsDeltaRho0(signalSimuArray,1,"Double_Non_Mon_");
