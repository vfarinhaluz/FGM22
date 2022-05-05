%This gives us the pdf of the lognormal distribution that is based on Einav
%et al (2013) paper.

function pdfValue = pdfEinav(mu, rho, rho_0)

    %Covariance matrix of logarithms (from Einav's paper):
    sigmaMu2=0.2;
    sigmaMuRho=-0.12;
    sigmaRho2=0.25;
    epsilon=10^-5;

%     The exact log averages from their paper are:
%     E_MuLog= .8; %8.2756;
%     E_RhoLog= 2; %-6.3909;

%   Fine-tuned averages:
    E_MuLog= 1.5;
    %When changing the grid for risk aversion, we adequately adjust the
    %distribution mean:
    if nargin==3
        E_RhoLog = log(rho_0) - .5*sigmaRho2;
    else
        E_RhoLog= 2;
    end
    
    %Calculatind pdf from local variation in CDF:
    pdfValue =(  mvncdf([log(mu + epsilon) ; log(rho + epsilon) ], [E_MuLog; E_RhoLog], [sigmaMu2 sigmaMuRho; sigmaMuRho sigmaRho2] ) -...
        mvncdf( [log(mu); log(rho)], [E_MuLog; E_RhoLog], [sigmaMu2 sigmaMuRho; sigmaMuRho sigmaRho2] )  ) / epsilon^2;
end


% The values for the grid should be:

% mu_0 = exp( E_MuLog + 0.5*sigmaMu2);
% muH: mu_0* exp ( (sigmaMu2)^.5 );
% muL: mu_0* exp ( - (sigmaMu2)^.5 );
% 
% rho_0 = exp( E_RhoLog + 0.5*sigmaRho2);
% rho_H: rho_0* exp ( (sigmaRho2)^.5 );
% rho_L: rho_0* exp ( - (sigmaRho2)^.5 );
