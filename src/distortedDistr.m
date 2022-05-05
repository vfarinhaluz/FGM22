%{
This code constructs a distribution for the model with binary risk aversion
levels. The distribution has a certain pattern of correlation designed to
induce negative welfare effects of a non-monotonic signal. 
Based on the analytical characterization of the paper, this distribution is
such that the value of w_0 * (1-w_0) is increasing in risk (see the paper
for more details on what the function w_0 represents).
%}
function pdfValue = distortedDistr(mu, rho, muL, muH, rhoL, rhoH)

    omega = 0.5 * (mu - muL) / (muH - muL);
    
    if rho==rhoL
        pdfValue = omega / (muH - muL);
    elseif rho==rhoH
        pdfValue = ( 1 - omega ) / (muH - muL);
    else
        error("This distribution only works for binary-continuous case.")
    end

end
