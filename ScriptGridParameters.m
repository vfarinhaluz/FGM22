% Code to go over grid of parameters and check signal effect.
addpath('./src');

clear;close all; clc;

%% Baseline parameters:


nContracts = 250
nRisk = 500
nRiskAv = 2

muL = 2;
muH = 4;
% rho_0=1;
% delta=1;
% rhoL = rho_0-delta/2;
% rhoH = rho_0+delta/2;

distributionName='Uniform';
distribution_pdf = @(m,y) ones(size(m))./( (muH-muL) );

BehavMass = 0;
tolerance = 1e-4;
maxIterations = 1e4;
alpha = 0.001; % alpha defines the size of the step
hyperParameters = HyperParameters;
hyperParameters.Tolerance = tolerance;
hyperParameters.MaxIterations = maxIterations;
hyperParameters.Alpha = alpha;

% MONOTONIC SIGNAL:
signal_A_Fun=@(m) (m-muL)./(muH-muL); %Must be a vectorized function in [0,1]
% NON MONOTONIC SIGNAL
% signal_A_Fun=@(m) SignalSimulation.nonMonotSignal(muL,muH,m); %Must be a vectorized function in [0,1]

%Creating grids
grid_Delta_percent=linspace(.2,.9,40);
grid_rho0= linspace(0.5,8,40);

%Create array varying all options. 
% Array dimensions are:
% lines: rho_0 values
% second (column): delta values
[Array_Delta_percent, Array_rho0] = meshgrid(grid_Delta_percent, grid_rho0);

%% Running grid
signalSimuArray=SignalSimulation.empty; %Empty array to store all simulations
length_grid=numel(Array_rho0);
i=1;
while i<=length_grid
    rho_0=Array_rho0(i);
    delta=Array_Delta_percent(i).*rho_0;
    rhoL = rho_0-delta/2;
    rhoH = rho_0+delta/2;

%     % Normal distribution (.9 coefficient of correlation)
%     distributionName='Normal - coef. of correlation 0.9';
%     sigma = [4,-.9*delta*2;
%             -.9*delta*2,delta^2]; %Variance covariance
%     mean = [(muL+muH)/2, rho_0];
%     distribution_pdf = @(mu,rho) mvnpdf([mu rho],mean,sigma);

    model=Model(muL, muH, rhoL, rhoH, nRisk, nRiskAv, nContracts, distribution_pdf, distributionName, BehavMass);
    solver = Solver(model,hyperParameters);
    signalSimuArray(i)= simulateSignal(signal_A_Fun,solver,"min");
    sprintf('#######  PERCENT COMPLETE: %.5g  ######',100*i/length_grid)
    i=i+1;
end

saveNoHistory(signalSimuArray, 'simuArray' );

%% Creating a table object with all relevant information:
 
plotHeatMaps(signalSimuArray)
