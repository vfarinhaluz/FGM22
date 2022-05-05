%{
Project: Risk Classification in Insurance Markets with Risk and Preference
Heterogeneity, by Vitor Farinha Luz, Humberto Moreira, Piero Gottardi
Matlab research assistance: Pedro Melgare
April 2nd, 2022

This code computes an approximate competitive equilibrium for the two-dimensional
model discussed in the main paper, allowing for any finite number of contracts and 
possible type realizations in each dimension. 
%}
%% Set up Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

%NEEDED line to find source code in subfolder:
addpath(strcat(pwd,'/src'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL PARAMETERS
%%% Grid parameters

nContracts = 250;
nRisk = 500;
nRiskAv = 2;

%%% Type distribution parameters
% Midpoint values for range of types:
rho_0= 8.8421;
mu_0 = 5;
% Sets range of types in each dimension:
delta = 0.2 * rho_0;
rhoL = rho_0  - delta/2;
rhoH = rho_0  + delta/2;
delta_mu = 0.8 * mu_0;
muH = mu_0  + delta_mu/2;
muL = mu_0 - delta_mu/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Type distribution

% Type distribution - Einav et al (2013) - Lognormal
distributionName='Einav et al (2013)';
distribution_pdf = @(mu, rho) pdfEinavFlex(mu, rho, mu_0, rho_0);

% %Uniform:
% distributionName='Uniform';
% distribution_pdf = @(m,y) ones(size(m))./( (muH-muL)*(rhoH-rhoL) );

% Simulation Parameters
BehavMass = 0; % Mass of behavioral types
tolerance = 1e-4;
maxIterations = 1e4;
alpha = 0.001; % alpha defines the size of the step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Not necessary to change from here forward
% Constructing HyperParameters object (details in HyperParameters.m):

hyperParameters = HyperParameters;
hyperParameters.Tolerance = tolerance;
hyperParameters.MaxIterations = maxIterations;
hyperParameters.Alpha = alpha;

%Constructing Model object (details in Model.m):

model=Model(muL,muH,rhoL,rhoH, nRisk, nRiskAv, nContracts,distribution_pdf, distributionName,BehavMass);

%Constructing Solver object (details in Solver.m):

solver = Solver(model,hyperParameters);
% solver.model = model;
% solver.hyperParameters = hyperParameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solving the model:
simResult = findEq(solver,'min');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting:
plotPriceTypeAssignmentExtreme(simResult,1); 
plotPoolMass(simResult,1);
plotDemand(simResult,1);
plotExtremesRhoDemand(simResult,1);
plotExtremesRhoDemandPrice(simResult,1);
plotDemand(simResult,1);
plotDemandLevelCurves(simResult,1)
