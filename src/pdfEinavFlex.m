%This gives us the pdf of the lognormal distribution that is based on Einav
%et al (2013) paper. It allows for flexibility in adjusting the first
%moments of the distribution of logarithms of types.

function pdfValue = pdfEinavFlex(mu, rho, mu_0, rho_0)

    %Covariance matrix of logarithms (from Einav's paper):
    sigmaMu2=0.2;
    sigmaMuRho=-0.12;
    sigmaRho2=0.25;

    %Set epsilon to compute numerical derivatives
    epsilon=10^-5;

%     The exact log averages from their paper are:
%     E_MuLog= 8.2756;
%     E_RhoLog= -6.3909;

%   Fine-tuned averages:
    %When changing the grid for risk aversion, we adequately adjust the
    %distribution mean:
    if nargin==4
        E_RhoLog = log(rho_0) - .5*sigmaRho2;
        E_MuLog = log(mu_0) - .5*sigmaMu2;
    else
        E_RhoLog= 2; % Default values
        E_MuLog= 1.5;
        disp("Reverting to default average values in pdf calculations.")
    end
    
    %Moments:
    MeanLog = [E_MuLog; E_RhoLog];
    CovMatLog = [sigmaMu2 sigmaMuRho;...
                sigmaMuRho sigmaRho2];

    %Calculatind pdf from local second derivative of CDF:

    diffMuEpsilon = ( mvncdf([log(mu + epsilon) ; log(rho + epsilon) ], MeanLog, CovMatLog ) -...
        mvncdf( [log(mu) ; log(rho + epsilon)], MeanLog, CovMatLog )  ) / epsilon;
    
    diffMuZero = ( mvncdf([log(mu + epsilon) ; log(rho) ], MeanLog, CovMatLog ) -...
        mvncdf( [log(mu); log(rho)], MeanLog, CovMatLog )  ) / epsilon;

    diff2 = ( diffMuEpsilon - diffMuZero ) /epsilon;

    pdfValue = diff2;
end
