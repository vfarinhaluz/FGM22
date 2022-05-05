%{
    This file defines the class Solver and implements a simple iterative
    fixed point algorithm to search for fixed points of EquilibriumOperator
    (defined in Model.m), which by construction must be equilibria of our
    model
%}
 % Defines class Solver that is used to construct a price iteration
 % algorithm to find an equilibrium

classdef Solver
    properties 
        model Model
        hyperParameters HyperParameters
    end
    methods
        
        %Constructor:
        function solver=Solver(model,hyperParameters)
        solver.model = model;
        solver.hyperParameters = hyperParameters;
        end

        % This is the algorithm based on operator EquilibriumOperator,
        % implemented in Model.m. This algorithm is vectorized for
        % efficiency. The algorithm is a simple fixed point iteration of 
        % EquilibriumOperator.
        function simResult = findEq(solver, method, verbose)           
            arguments %Validating to make sure arguments work
                solver Solver
                method ="min" %last or min allowed
                verbose = 0;
            end
                
            timebegin=tic;
    
            tol = solver.hyperParameters.Tolerance;
            maxIt = solver.hyperParameters.MaxIterations;
            alpha = solver.hyperParameters.Alpha;
            modelTemp = solver.model;
            nRisk = length(modelTemp.RiskSet);
            nRiskAv = length(modelTemp.RiskAversionSet);
    
            %Initial guess of price function:
            rho0=(modelTemp.RiskAversionSet(nRiskAv)+modelTemp.RiskAversionSet(1))/2;
            muH = modelTemp.RiskSet(nRisk);
            pricefunction = @(x) x.*(rho0.*(1-x+log(x))+muH);
            p = pricefunction(modelTemp.ContractSpace);
        
            %Objects for iterations:
            errorHist = nan(maxIt,1);
            pHist=nan(length(modelTemp.ContractSpace), maxIt);
            nIt = 1;
            error = 1;
            errorMin=1;
            pMinError=p;
            msgWind = floor(maxIt/4);            
            tloop=tic;
            tloopel = 0;
            
            % Set maximum time to search for an equilibrium
            tloopmax = 1000;
    
            while (nIt <= maxIt && error > tol && tloopel <= tloopmax)
                tloopel = toc(tloop);
                tloop=tic;
                if verbose==1 && mod(nIt,msgWind)==0                    
                    error
                    nIt                    
                end
                [error,deltaPrice]=EquilibriumOperator(modelTemp,p);
                pHist(:,nIt)=p;
                errorHist(nIt) = error;
                if error<errorMin;
                    pMinError=p;
                    errorMin=error;
                end
                p = p + alpha.*deltaPrice;
                nIt = nIt + 1;
            end
            
            if tloopel > tloopmax
                disp('Error: maximum iteration time reached, no equilibrium found')
            end
    
            str=strcat('Total seconds to complete:  ',  num2str(toc(timebegin)), ' seconds');
            disp(str);
    
            %Saving output:
            simResult=SimulationResult;
            simResult.model=modelTemp;
            if method=="min";
                simResult.eqPrice=pMinError;
            elseif method=="last";
                simResult.eqPrice=p;
            else
                error('Incorrect method defined. Use either last or min!')
            end
            simResult.eqUtilityMatrix=CARA_UtilityMatrix(modelTemp,p);
            simResult.eqDemandMatrix=DemandMatrix(modelTemp,p);
            simResult.errorHist=errorHist;
            simResult.errorMin=errorMin;
            simResult.hyperParameters=solver.hyperParameters;

        end       
        
    end
end 