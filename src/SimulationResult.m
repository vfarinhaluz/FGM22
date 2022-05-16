%{
This file defines SimulationResult class, which collects equilibrium
objects. It also implements functions to plot equilibrium objects, such as 
prices, demand, and type assignment.
%}

classdef SimulationResult
    %This class represents the result of a simulation
    %   It stores all the simulated outcomes (price, utility, errors, etc)
    %   as well as the original model details from the simulation.
    
    properties
        model Model
        eqPrice (:,1) double
        eqUtilityMatrix double
        eqDemandMatrix double
        errorHist double
        hyperParameters HyperParameters
        errorMin double
    end
    
    methods
        % This function plots a error x iteration graph
        function plotErrorHist(simResult,save_indicator)
            figure;
            plot(simResult.errorHist);
            title('Error');
            xlabel('Iteration');
            ylabel('Error');
            if nargin == 2 & save_indicator==1
                saveas(gcf,'figures/errorHist.png')
            end
        end

        % This function plots the equilibrium price function
        function plotEquilPrice(simResult,save_indicator)
            figure
            plot(simResult.model.ContractSpace,simResult.eqPrice)
            title('Equilibrium price function');
            xlabel('x')
            ylabel('p')
            axis([0 1 -Inf inf]);
            if nargin == 2 && save_indicator==1
                saveas(gcf,'figures/equilibriumPriceFunction.png')
            end
        end

        % This function compares 1. the equilibrium price function found by
        % our iteration algorithm for a model with heterogeneity in both
        % risk aversion and risk level and 2. the equilibrium price
        % function for a model with one dimensional heterogeneity 
        % on the risk level (in this case, there is a simple analytical
        % solution)
        function plotCompareUniDimension(simResult,save_indicator)
            muL=simResult.model.muL;
            muH=simResult.model.muH;
            rhoL=simResult.model.rhoL;
            rhoH=simResult.model.rhoH;
            contracts=simResult.model.ContractSpace;
    
            rho0 = (rhoH-rhoL)/2;
            pUniFun = @(x) x.*(rho0.*(1-x+log(x))+muH);
            pUniGrid = pUniFun(simResult.model.ContractSpace);
    
            figure;
            plot(contracts,pUniGrid,'red','DisplayName','Analytical');
            hold on;
            plot(contracts,simResult.eqPrice,'blue','DisplayName','Numerical');
            title('Two dimensional vs one dimensional solution');
            legend('One dimensional','Two dimensional');
            hold off;

            if nargin == 2 && save_indicator==1
                saveas(gcf,'figures/Compare_1dim.png')
            end
        end

        % This function compares the average cost with the price for each
        % contract. In an equilibrium, the two should exactly match.
        function plotCompareAverageCost(simResult)
            [AC,~,~] = AverageCost(simResult.model,simResult.eqPrice);
            contracts = simResult.model.ContractSpace;
            figure;
            plot(contracts,simResult.eqPrice,'DisplayName','Price');
            hold on;
            plot(contracts,AC,'DisplayName','Average Cost');
            title('Price and Average Cost');
            legend('Price','Average Cost');
            hold off;
            saveas(gcf,'figures/PriceAverageCost.png')
        end

        % This function plots the demand for consumers with the lowest and
        % highest levels of risk aversion.
        function plotExtremesRhoDemand(simResult,save_indicator)
            muL=simResult.model.muL;
            muH=simResult.model.muH;
            nRiskAv=length(simResult.model.RiskAversionSet);
            demandMat=DemandMatrix(simResult.model,simResult.eqPrice);
            
            figure;
            riskSetGraph = linspace(muL,muH,1000);
            DL = interp1(simResult.model.RiskSet,demandMat(:,1),riskSetGraph);
            DH = interp1(simResult.model.RiskSet,demandMat(:,nRiskAv),riskSetGraph);
            plot(riskSetGraph,DL);
            hold on;
            plot(riskSetGraph,DH);
            legend('rhoL','rhoH');
            xlabel('Risk level');
            ylabel('Coverage level');
            title('Type allocation for rhoL, rhoH');
            hold off;

            if nargin == 2 && save_indicator==1
            saveas(gcf,'figures/extremesRhoDemand.png')
            end
            
        end

        % This function does the same as plotExtremesRhoDemand, but also
        % includes the price function in the figure
        function plotExtremesRhoDemandPrice(simResult,save_indicator)
            plotExtremesRhoDemand(simResult);
            hold on;
            [~,isTraded]= AverageCostUpdate(simResult.model, simResult.eqPrice);
            plot(simResult.eqPrice(isTraded) ./ simResult.model.ContractSpace(isTraded)...
                , simResult.model.ContractSpace(isTraded))
            legend('rhoL','rhoH','price');
        
            if nargin == 2 & save_indicator==1
                saveas(gcf,'figures/extremesRhoDemand.png')
            end    
        end

        % This function plots the type assignment function for consumers
        % with the highest and lowest risk aversion values
        function plotPriceTypeAssignmentExtreme(simResult,save_indicator, isContinuous)

            arguments
                simResult
                save_indicator =0
                isContinuous = 0
            end
            model0=simResult.model;
            price=simResult.eqPrice;
            contracts=simResult.model.ContractSpace;
            riskSet=simResult.model.RiskSet;
            
            demandMat=DemandMatrix(model0, price);
            [~,isTraded]= AverageCostUpdate(model0, price);
            
            lineWidth=2;
            font=15;
            fig=figure;
            hold on;
            plot(demandMat(:,1), riskSet,'Color','red', 'LineWidth',lineWidth);
            plot(demandMat(:,end), riskSet,'Color','blue', 'LineWidth',lineWidth);
            plot(contracts(isTraded), price(isTraded) ./ contracts(isTraded), 'Color','black', 'LineWidth',lineWidth);
            xlabel('Coverage level');
            ylabel('Risk level');
            if isContinuous==1
                legend("Type assignment – m(\cdot," + round(model0.rhoL,1) + ")","Type assignment – m(\cdot," + round(model0.rhoH,1) + ")",'Location','southeast','FontSize', font);
                title('Type assignments and price per unit', FontSize=font);
            else
                legend('Type assignment – m_l','Type assignment – m_h','Price per unit','Location','southeast','FontSize', font);
                title('Type assignments and price per unit', FontSize=font);
            end
            set(gca,'Fontsize', 16);
            hold off;

            if save_indicator==1
            f=fullfile('./figures','SimuTypeAssign.pdf');
            saveas(fig,f);
            end

        end

        % This function makes a 3D plot displaying the coverage demanded
        % by each consumer type
        function plotDemand(simResult,save_indicator)
            demandMat=DemandMatrix(simResult.model,simResult.eqPrice);
        
            fig=figure;
            surf(simResult.model.RiskAversionSet,simResult.model.RiskSet,demandMat);
            ylabel('Risk level');
            xlabel('Risk aversion');
            zlabel('Coverage level');
            title('Demand levels');

            if nargin == 2 & save_indicator==1
                saveas(fig,'figures/Allocation3D.png');
            end            
        end
    
        % This function plots level curves for the demand function on the
        % risk aversion x risk level plane
        function plotDemandLevelCurves(simResult,save_indicator)   
            arguments
                simResult
                save_indicator = 0;
            end
            
            fig=figure; hold on;
            x = simResult.model.RiskAversionSet;
            y = simResult.model.RiskSet;
            [X,Y] = meshgrid(x,y);
            Z = DemandMatrix(simResult.model,simResult.eqPrice);
            contour(X,Y,Z,'ShowText','on');
            ylabel('Risk level');
            xlabel('Risk aversion');
            title('Level curves of demand function');

            if save_indicator==1
                saveas(fig,'figures/DemandLevelCurves.png');
            end
        end

        % This function creates a heatmap displaying the demand for each
        % consumer type
        function plotDemandHeatmap(simResult,save_indicator)   
            arguments
                simResult
                save_indicator = 0;
            end
        
            vec_ra = simResult.model.RiskAversionSet;
            vec_r = simResult.model.RiskSet;
            [ra_grid, risk_grid] = meshgrid(vec_ra,vec_r);
            demand_grid = DemandMatrix(simResult.model,simResult.eqPrice);    
    
            tableDemand=table;
            tableDemand.riskAversion=round(ra_grid(:),1);
            tableDemand.risk=round(risk_grid(:),1);
            tableDemand.demand=round(demand_grid(:), 2);
    
            figure;
            p=heatmap(tableDemand,"riskAversion","risk","ColorVariable","demand","Colormap",jet,...
                CellLabelColor = 'none', FontSize=14);
            riskaversion1 = tableDemand.riskAversion(1); %This is a hack to reorder heatmap appropriately.
            sorty(p,string(riskaversion1), 'descend');
            title('Equilibrium demand heatmap');
            xlabel('Risk averion \rho');
            ylabel('Risk level \mu');

            if save_indicator==1
                saveas(p,'figures/DemandHeatmap.pdf');
            end
        end

        % This function makes a 3D plot showing the utility level for each
        % risk aversion and risk level pair
        function plotUtility3D(simResult,save_indicator)    
            fig=figure;
            x = simResult.model.RiskAversionSet;
            y = simResult.model.RiskSet;
            [X,Y] = meshgrid(x,y);
            Z = CARA_UtilityMatrix(simResult.model,simResult.eqPrice);
            surf(X,Y,Z);
            ylabel('Risk level');
            xlabel('Risk aversion');
            title('Equilibrium Utility levels');

            if save_indicator==1
                saveas(fig,'figures/Utility3D.png');
            end
        end

        % This function plots the mass of consumers choosing each contract
        % in equilibrium.
        function plotPoolMass(simResult,save_indicator)   
            arguments
                simResult
                save_indicator =0;
            end
            
            Contracts = simResult.model.ContractSpace;
            massInEachContract = poolMass(simResult.model, simResult.eqPrice);
            fig=figure; hold on;
            bar(Contracts, massInEachContract);
            ylabel('Mass of agents');
            xlabel('Coverage');
            set(gca,'Fontsize', 16);
            if length(Contracts)<=30; xticks(Contracts); end
            xtickformat('%.2f');
            title('Mass of agents in each pool');
        
            if save_indicator==1
                saveas(fig,'figures/massPool.pdf');
            end
        end

        % tableGridEqSimulation stores some variables of interest in a
        % table for a grid of simulations
        function simuTable = tableGridEqSimulation(eqSimuArray)
            for i=1:length(eqSimuArray)
                vectorContractNumber(i)=length(eqSimuArray(i).model.ContractSpace);
                vectorRiskNumber(i)=length(eqSimuArray(i).model.RiskSet);
                Array_errors(i)= eqSimuArray(i).errorMin;
            end
                simuTable=table;
                simuTable.row=(1:length(vectorContractNumber))';
                simuTable.contractNumber=vectorContractNumber(:);
                simuTable.riskNumber=vectorRiskNumber(:);
                simuTable.errors=Array_errors(:);
        end        

    end
end
