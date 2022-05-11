%{
    This file defines the class SignalSimulation, solves a binary signal
    disclousure exercise, compute welfare effects, plots changes in
    equilibrium objects such as prices and utility levels.
%}

classdef SignalSimulation
%     This class saves the results of a simulation to study the effect of
%     the signal release. Hence it contains three equilibria with different
%     distributions: prior, posterior A, and posterior B. It also contains
%     the function handle for the signalA probability function.

    properties 
        signalFun
        prior SimulationResult
        signalA SimulationResult
        signalB SimulationResult
    end
    
    methods
        % SignalSimulation constructs an object of the class
        % SignalSimulation
        function sim=SignalSimulation(simulation_prior,simulation_A,simulation_B,signal_A_fun)
            sim.prior=simulation_prior;
            sim.signalA=simulation_A;
            sim.signalB=simulation_B;
            sim.signalFun=signal_A_fun;
        end
	          
        % UtilCompareSignal does some welfare computations on the
        % disclousure of a signal. It outputs the average change for each
        % consumer type (averaged out of all possible signal realizations),
        % the share of the population that gains from the disclousure of
        % the signal and the aggregate welfare impact (with uniform welfare
        % weights) of the disclousure
        function [UtilChangeAverageMat, sharePopulationGain, AggregWelfareImpact]=UtilCompareSignal(signalSimulation,ignoreTooSmall)            
            arguments
                signalSimulation SignalSimulation
                ignoreTooSmall double = 0
            end
            
            riskSet=signalSimulation.prior.model.RiskSet;
            RiskAvSet=signalSimulation.prior.model.RiskAversionSet;
            nRiskAvSet=length( RiskAvSet );
            Distrib=signalSimulation.prior.model.Distribution;
            signal_A_Fun=signalSimulation.signalFun;
            
            %Prices:
            price_0=signalSimulation.prior.eqPrice;
            price_A=signalSimulation.signalA.eqPrice;
            price_B=signalSimulation.signalB.eqPrice;
            
            %Utility matrices:
            UtilMatrixPrior=CARA_UtilityMatrix(signalSimulation.prior.model,price_0);
            UtilMatrixSignalA=CARA_UtilityMatrix(signalSimulation.signalA.model,price_A);
            UtilMatrixSignalB=CARA_UtilityMatrix(signalSimulation.signalB.model,price_B);
            
            %Average utility:
            probSignalA_Matrix=repmat(signal_A_Fun(riskSet),1,nRiskAvSet);
            UtilMatrixAverage=UtilMatrixSignalA.*probSignalA_Matrix +...
                        UtilMatrixSignalB .* (1 - probSignalA_Matrix);
            UtilChangeAverageMat=UtilMatrixAverage-UtilMatrixPrior;

            sharePopulationGain= sum(...
                (UtilChangeAverageMat>=0).* Distrib,"all");

            %Aggregate welfare impact:
            AggregWelfareImpact=sum(...
                UtilChangeAverageMat .* Distrib...
                ,'all');

            if ignoreTooSmall==1 %If we want to eliminate changes that are too small
                error=signalSimulation.prior.errorHist;
                error=error(~isnan(error));
                thresholdIgnore=error(end);
                isvalid=(  abs(UtilChangeAverageMat)  >  thresholdIgnore  );
            
                UtilChangeAverageMat=UtilChangeAverageMat.* isvalid;

                sharePopulationGain= sum(...
                (UtilChangeAverageMat>=0).* Distrib,"all");

            end            
        end 
     

        % tableGridSimulation generates several tables with values of
        % relevant variables for a grid of simulations, including: 1. the
        % errors in the iteration algorithm, 2. the share of the consumer
        % population which receives an interim utility increase after the signal
        % disclousure, 3. the aggregate welfare effect of the signal
        % disclousure, 4. the average equivalent variation in the
        % population, 5. the variation in risk and risk aversion in the
        % population, 6. the mid-value for risk aversion.
        function simuTable = tableGridSimulation(signalSimuArray)
            for i=1:length(signalSimuArray)
                Array_rho0(i)=mean(signalSimuArray(i).prior.model.RiskAversionSet);
                Array_rho0(i)= round(Array_rho0(i), 10);
                Array_Delta_percent(i)=...
                    range( signalSimuArray(i).prior.model.RiskAversionSet) ./Array_rho0(i);
                Array_Delta_percent(i)= round(Array_Delta_percent(i), 5);
                Array_errors(i)= max( [signalSimuArray(i).prior.errorMin ...
                                        signalSimuArray(i).signalA.errorMin ...
                                        signalSimuArray(i).signalB.errorMin]  );
                [~,Array_sharePopSignal(i),Array_AggregateWelfSignal(i)]= UtilCompareSignal(signalSimuArray(i));
                [~, Array_AverageEquivVariation(i)] = signalEquivalentVariation(signalSimuArray(i));
                Array_MidRisk(i)=mean(signalSimuArray(i).prior.model.RiskSet);
                Array_deltaRiskPercent(i)=range(signalSimuArray(i).prior.model.RiskSet) ./ Array_MidRisk(i);
                Array_deltaRiskPercent(i)= round(Array_deltaRiskPercent(i), 5);
            end
            simuTable=table;
            simuTable.row=(1:length(Array_rho0))';
            simuTable.delta_percent=round(Array_Delta_percent(:), 2);
            simuTable.rho0=round(Array_rho0(:), 2);
            simuTable.errors=Array_errors(:);
            simuTable.shareWelfare=round(Array_sharePopSignal(:), 1);
            simuTable.AggregateWelfSignal=Array_AggregateWelfSignal(:);
            simuTable.AverageEquivVariation = Array_AverageEquivVariation(:);
            simuTable.deltaRisk_percent=Array_deltaRiskPercent(:);
        end

        % This function generates three heatmaps plots for a grid of delta (risk
        % aversion dispersion) and rho0 (risk aversion mean). The plots
        % display 1. the error term in the iteration algorithm, 2. the
        % mass of consumers with interim improvement, and 3. the aggregate
        % welfare effect of the signal disclousure
        function [p1, p2, p3] = plotHeatMapsDeltaRho0(signalSimulationArray, save_indicator, titleDescription)            

            arguments
                signalSimulationArray SignalSimulation
                save_indicator = 0
                titleDescription = ""
            end
            
            results = tableGridSimulation(signalSimulationArray);
        
            p1=figure;
            heatmap(results,"rho0","delta_percent","ColorVariable","errors","Colormap",flipud(autumn),...
                CellLabelFormat = '%.2E',FontSize=14);
            title('Final error');

            figure;
            p2=heatmap(results,"rho0","delta_percent","ColorVariable","shareWelfare","Colormap",autumn,...
                CellLabelFormat = '%.1g',FontSize=14);
            title("Mass with interim improvement " + titleDescription);
            xlabel('\rho_{0}');
            
            ylabel('\delta/\rho_{0}'); %,'Interpreter','latex')
            figure;
            p3=heatmap(results,"rho0","delta_percent","ColorVariable","AggregateWelfSignal","Colormap",autumn,...
                CellLabelFormat = '%.2E',FontSize=14,Title="Aggregate welfare gain from signal");
            title('Signal aggregate welfare effect');
            
           if save_indicator
            saveas(p1,'figures/HeatmapErrors.pdf');
            saveas(p2,'figures/HeatmapWelfare.pdf');
            saveas(p3,'figures/HeatmapAggWelfare.pdf');
           end        
        end
    
        % This function does the same than plotHeatMapsDeltaRho0, but for a
        % grid of risk level dispersion and rho0 (mean risk aversion)
        function plotHeatMapsRiskDelta(signalSimulationArray, save_indicator)            
            arguments
                signalSimulationArray SignalSimulation
                save_indicator = 0;
            end 
                   
            results = tableGridSimulation(signalSimulationArray);
        
            p1=figure;
            heatmap(results,"deltaRisk_percent","delta_percent","ColorVariable","errors","Colormap",flipud(autumn),...
                CellLabelFormat = '%.2E',FontSize=14);
            title('Final error');
            ylabel('\delta/\rho_{0}');
            xlabel('\Delta\mu/\mu_{0}');
            figure;
            p2=heatmap(results,"deltaRisk_percent","delta_percent","ColorVariable","shareWelfare","Colormap",autumn,...
                CellLabelFormat = '%.1g',FontSize=14);
            title('Mass with interim improvement');
            ylabel('\delta/\rho_{0}');
            xlabel('\Delta\mu/\mu_{0}');
            figure;
            p3=heatmap(results,"deltaRisk_percent","delta_percent","ColorVariable","AggregateWelfSignal","Colormap",autumn,...
                CellLabelFormat = '%.2E',FontSize=14,Title="Aggregate welfare gain from signal");
            title('Signal aggregate welfare effect');
            ylabel('\delta/\rho_{0}');
            xlabel('\Delta\mu/\mu_{0}');
            
           if save_indicator
            saveas(p1,'figures/HeatmapErrors.pdf');
            saveas(p2,'figures/HeatmapWelfare.pdf');
            saveas(p3,'figures/HeatmapAggWelfare.pdf');
           end        
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PLOTS
        % This function produces several plots comparing the equilibrium
        % under the prior distribution with the equilibria after the signal
        % disclousure. It compares prices and welfare. It also plots the
        % signal function.
        function plotSignalSimulation(signalSimulation,save_indicator)
            arguments
                signalSimulation SignalSimulation
                save_indicator = 0;                
            end

            % Load model parameters
            muL=signalSimulation.prior.model.muL;
            muH=signalSimulation.prior.model.muH;
            rhoL=signalSimulation.prior.model.rhoL;
            rhoH=signalSimulation.prior.model.rhoH;        
            riskSet=signalSimulation.prior.model.RiskSet;
            RiskAvSet=signalSimulation.prior.model.RiskAversionSet;
            contractSpace=signalSimulation.prior.model.ContractSpace;
            nRiskAvSet=length(RiskAvSet);
            signal_A_Fun=signalSimulation.signalFun;
            
            %Prices:
            price_0=signalSimulation.prior.eqPrice;
            price_A=signalSimulation.signalA.eqPrice;
            price_B=signalSimulation.signalB.eqPrice;
            
            %Prices per unit:
            Unit_price_0=signalSimulation.prior.eqPrice./contractSpace;
            Unit_price_A=signalSimulation.signalA.eqPrice./contractSpace;
            Unit_price_B=signalSimulation.signalB.eqPrice./contractSpace;
            
            % Utility matrices:
            UtilMatrixPrior=CARA_UtilityMatrix(signalSimulation.prior.model,price_0);
            UtilMatrixSignalA=CARA_UtilityMatrix(signalSimulation.signalA.model,price_A);
            UtilMatrixSignalB=CARA_UtilityMatrix(signalSimulation.signalB.model,price_B);

            %Average utility:
            probSignalA_Matrix=repmat(signal_A_Fun(riskSet),1,nRiskAvSet);
            UtilMatrixAverage=UtilMatrixSignalA.*probSignalA_Matrix +...
                        UtilMatrixSignalB .* (1 - probSignalA_Matrix);
            UtilChangeAverage=UtilMatrixAverage-UtilMatrixPrior;

            demandMat = DemandMatrix(signalSimulation.prior.model,Unit_price_0.*contractSpace);
            xmin=min(demandMat,...
                [],"all")-.01;

            %% Plotting
            % Prices            
            pricesPlot = figure;
            plot(contractSpace,price_0,'LineWidth',2); 
            hold on;
            plot(contractSpace,price_A);
            plot(contractSpace,price_B);    
            fontSize=14;
            title('Price comparison - Prior and both signals',...
                'FontSize',fontSize);
            axis([0 1 -Inf inf]);
            xlabel('coverage','FontSize',fontSize);
            ylabel('prices','FontSize',fontSize);
            legend('Prior','Signal A', 'Signal B',...
                'FontSize',fontSize,...
                'Location','southeast'); 
            hold off;

            % Price per unit      
            pricePerUnitPlot=figure; 
            hold on;
            set(pricePerUnitPlot,'defaultLineLineWidth',2);
            set(gca,'Fontsize', 16);
            plot(contractSpace,Unit_price_0); 
            plot(contractSpace,Unit_price_A);
            plot(contractSpace,Unit_price_B);    
            title('Price per unit comparison - Prior and both signals')
            axis([xmin 1 muL muH]);
            xlabel('Coverage (x)');
            ylabel('Prices per unit (p/x)');
            legend('No signal','Signal A', 'Signal B',...
                'Location','southeast'); 
            hold off;

            %Signal function
            signalFunctionPlot=figure;
            hold on;
            fplot(@(m) signal_A_Fun(m),[muL muH]);
            title('Probability of signal A, conditional on risk ( P(s=A|μ) )')
            ylabel('Probability of signal A')
            xlabel('Risk Level')
            set(gca,'Fontsize', 16, 'linewidth', 1)
            hold off;
            
            %Utility change
            uChangePlot=figure;
            mesh(RiskAvSet,riskSet,UtilChangeAverage);
            hold on;
            title('Average utility change from signal - CARA utility')
            ylabel('Risk level')
            xlabel('Risk aversion')
            zlabel('Utility change')
            set(gca,'Fontsize', 16, 'linewidth', 1)
            hold off;

            %Set of types with positive utility change:
            setUChanges = figure;
            RiskMatrix=repmat(riskSet,1,nRiskAvSet);
            RiskAvMatrix=repmat(RiskAvSet',length(riskSet),1);
            scatter( RiskMatrix(UtilChangeAverage>=0), ...
                RiskAvMatrix(UtilChangeAverage>=0),...
                'filled' );
            axis([muL muH rhoL rhoH]);
            title('Set of types with interim improvement - CARA utility');
            xlabel('Risk level');
            ylabel('Risk aversion');
            set(gca,'Fontsize', 16, 'linewidth', 1);
            
            % Equivalent variation extremes
            [EV_Mat, ~] = signalEquivalentVariation(signalSimulation);
            
            plotEquivVariation=figure;
            hold on;
            set(gca,'Fontsize', 16);
            set(plotEquivVariation,'defaultLineLineWidth',2);
            riskSetGraph = linspace(muL,muH,1000);
            EVL = interp1(riskSet,EV_Mat (:,1), riskSetGraph);
            EVH = interp1(riskSet,EV_Mat (:,nRiskAvSet), riskSetGraph);
            plot(riskSetGraph,EVL);
            plot(riskSetGraph,EVH);
            xlim([0.9*muL 1.05*muH])
            hold off;
            legend('\rho_{L}','\rho_{H}');
            xlabel('Risk level (\mu)');
            ylabel('Equivalent variation');
            title('Equivalent variation of signal disclosure:');

            %% Saving plots
            if nargin == 2 && save_indicator==1
                 saveas(pricesPlot,'figures/plotPrices.pdf')
                 saveas(pricePerUnitPlot,'figures/plotPricesPerUnit.pdf')
                 saveas(signalFunctionPlot,'figures/signalFunctionPlot.pdf')
                 saveas(uChangePlot,'figures/uChangePlot.pdf')
                 saveas(setUChanges,'figures/setUChanges.pdf')
                 saveas(plotEquivVariation,'figures/equivalentVariation.pdf')
            end
        end            
               
        % signalEquivalentVariation computes how much each consumer would
        % be willing to pay for the signal to be disclosed (before knowing
        % its particular realization of the signal). It also computes the
        % mean willingness to pay across consumers
        function [EV, meanEV] = signalEquivalentVariation(signalSimulation)    
            arguments
                signalSimulation SignalSimulation
            end            
        
            riskSet=signalSimulation.prior.model.RiskSet;
            nriskSet = length(riskSet);
            riskAvSet = signalSimulation.prior.model.RiskAversionSet;
            riskAvMat = repmat(riskAvSet', nriskSet, 1);
            nRiskAvSet=length( riskAvSet );
            signal_A_Fun=signalSimulation.signalFun;
            distr = signalSimulation.prior.model.Distribution;
            
            %Prices:
            price_0=signalSimulation.prior.eqPrice;
            price_A=signalSimulation.signalA.eqPrice;
            price_B=signalSimulation.signalB.eqPrice;
            
            %Utility matrices:
            UtilMatrixPrior=CARA_UtilityMatrix(signalSimulation.prior.model,price_0);
            UtilMatrixSignalA=CARA_UtilityMatrix(signalSimulation.signalA.model,price_A);
            UtilMatrixSignalB=CARA_UtilityMatrix(signalSimulation.signalB.model,price_B);
            
            %Average utility:
            probSignalA_Matrix=repmat(signal_A_Fun(riskSet),1,nRiskAvSet);
            UtilMatrixAverage=UtilMatrixSignalA.*probSignalA_Matrix +...
                        UtilMatrixSignalB .* (1 - probSignalA_Matrix);
            
            CEsignal = - log( - UtilMatrixAverage ) ./ riskAvMat;
            CEzero = - log( - UtilMatrixPrior ) ./ riskAvMat;
            EV = CEsignal - CEzero;
            meanEV = sum(EV .* distr, 'all');
        end
        

    end

    methods (Static)

        % nonMonotSignal generates a non-monotonic signal function
        % This function uses the model parameters to construct a quadratic
        % non-monotonic signal. It has two optional parameters that determine how
        % informative the signal is.
        function probSignalA=nonMonotSignal(muL,muH,m,plow,phigh)      
            arguments
                muL double
                muH double
                m double   
                plow =.2
                phigh =.8      
            end        
                

            b=phigh;
            a=-(phigh-plow)*4/( (muH-muL)^2 );
            probSignalA = a.*(m-(muH+muL)/2).^2 +b;
            probSignalA = probSignalA .* (m>=muL) .*(m<=muH);
        end

        % This function throws away the error History from the signal simulation
        % array and saves it. The goal is to store a smaller object.
        function saveNoHistory(signalSimuArray,filename)            
            arguments
                signalSimuArray SignalSimulation
                filename char ='unamedFile'             
            end            
            
            for i=1:length(signalSimuArray)
                signalSimuArray(i).prior.errorHist=NaN;
                signalSimuArray(i).signalA.errorHist=NaN;
                signalSimuArray(i).signalB.errorHist=NaN;
            end
            save(filename,'signalSimuArray','-v7.3')
        end

         % The function updateModelSignal creates a new Model object with a
         % type distribution derived from a Bayesian update after observing
         % signal signal_Fun
        function modelNew = updateModelSignal(model,signal_Fun,NewDistrName)
            %Removing critical objects:
            muL=model.muL;
            muH=model.muH;
            rhoL=model.rhoL;
            rhoH=model.rhoH;
            nRisk=length(model.RiskSet);
            nRiskAv=length( model.RiskAversionSet );
            nContracts=length( model.ContractSpace );
            distribution_pdf= model.distribution_pdf;
            BehavMass=model.BehavioralMass;
	        %UtilHandle=model_temp.U;        
        
            %Signal A:
            %This new density is not properly normalized, but this is not an issue
            %since the Model constructor normalizes the distribution when
            %constructing the new disrtibution in the grid of types.
        
            distribution_pdf_Signal=@(m,r) distribution_pdf(m,r) .* signal_Fun(m);
            distributionName=NewDistrName;    
            modelNew=Model(muL,muH,rhoL,rhoH, nRisk, nRiskAv, nContracts,distribution_pdf_Signal, distributionName,BehavMass);
        end

        % simulateSignal solves a model with disclousure of a binary
        % signal, with probability of signal 'A' being received for a
        % risk type m given by signal_A_fun(m)
        function signalSimulation=simulateSignal(signal_A_Fun,solver, method)
            arguments
                signal_A_Fun {mustBeA(signal_A_Fun,'function_handle')}
                solver Solver
                method ="min"
            end    
            tstart=tic;    
    
            %Prior
            eq_prior=findEq(solver,method);
            
            %Signal A:
            %This new density is not properly normalized, but this is not an issue
            %since the Model constructor normalizes the distribution when
            %constructing the new disrtibution in the grid of types.
            solver_temp=solver;
            model_temp=solver_temp.model;
            model_temp=SignalSimulation.updateModelSignal(model_temp,signal_A_Fun,'Signal A');
            solver_temp.model=model_temp;
            eq_signalA=findEq(solver_temp,method); %Storing
        
            %Signal B:
            solver_temp=solver;
            model_temp=solver_temp.model;
            model_temp=SignalSimulation.updateModelSignal(model_temp,@(m) (1 - signal_A_Fun(m) ),'Signal B');
            solver_temp.model=model_temp;
            eq_signalB=findEq(solver_temp,method); %Storing
    
            signalSimulation=SignalSimulation(eq_prior,eq_signalA,eq_signalB,signal_A_Fun);
    
            str=strcat('Signal simulation total seconds to complete:  ',  num2str(toc(tstart)), ' seconds');
            disp(str);
        end
    end
end




