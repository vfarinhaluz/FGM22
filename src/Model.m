%{
This file defines the class Model, which collects model parameters and supporting objects, and functions to:
1. construct a Model object
2. compute an 3d array which stores utility levels for all consumer types and all contracts
3. compute a matrix which stores the maximum attainable utility level for each consumer type
4. compute demand for all consumers types
5. compute average cost for each contract
6. update a guess for the price vector (used in Solver.m findEq function)
7. compute the mass of consumers choosing each contract
%}
classdef Model
    % The class Model collects model parameters and supporting objects.
    properties 
        ContractSpace (:,1)
        RiskSet (:,1)
        RiskAversionSet (:,1)
        BehavioralMass (1,1) double
        Distribution double
        muL double
        muH double
        rhoL double
        rhoH double
        U_Gross_Array double 
        distributionName char
        distribution_pdf function_handle
    end
    
    methods
	
	    % Model constructs a Model object
        function model=Model(muL,muH,rhoL,rhoH, nRisk, nRiskAv, nContracts,distribution_pdf, distributionName,BehavioralMass)
            model.ContractSpace=linspace(.01,.99,nContracts)';
            model.RiskSet=linspace(muL,muH,nRisk)';
            model.RiskAversionSet=linspace(rhoL,rhoH,nRiskAv)';
            model.distributionName=distributionName;
            model.distribution_pdf=distribution_pdf;
            model.BehavioralMass=BehavioralMass;
            model.muL=muL;
            model.muH=muH;
            model.rhoL=rhoL;
            model.rhoH=rhoH;

            %Computing distribution in grid:
            sumTemp = 0;
            for i = 1:nRisk
                for j = 1:nRiskAv
                    phi(i,j) = distribution_pdf(model.RiskSet(i), model.RiskAversionSet(j));
                    sumTemp = sumTemp + phi(i,j);
                end
            end
            phi = phi./sumTemp;

            model.Distribution=phi;
            model=ComputeUtilArray(model);
        end

	    % ComputeUtilArray creates a 3-dim array with all utility levels for all contracts and types:
	    % Dimensions represent, in order: price, risk and risk aversion
        function modelNew=ComputeUtilArray(model)
            modelNew=model;
            nRisk = length(model.RiskSet);
            nRiskAv = length(model.RiskAversionSet);
            nContracts = length(model.ContractSpace);            
	             
            contractArray=repmat(model.ContractSpace,...
                1,nRisk,nRiskAv);
            riskArray = repmat(model.RiskSet',nContracts,1,nRiskAv);
            riskAvArray= repmat( reshape(model.RiskAversionSet,1,1,nRiskAv)...
                ,nContracts,nRisk,1); 
            uGrossArray= riskArray.*contractArray -(riskAvArray./2).*( (1-contractArray).^2);
            
            modelNew.U_Gross_Array=uGrossArray;
        end

        % CARA_UtilityMatrix computes maximum utility for each consumer type given a price vector p
        function utilMat=CARA_UtilityMatrix(model,p) %Correct one to use!
            nRisk = length(model.RiskSet);
            nRiskAv = length(model.RiskAversionSet);
            nContracts = length(model.ContractSpace);
            uGrossArray=model.U_Gross_Array;

            RiskAvArray_2D = repmat( model.RiskAversionSet', nRisk,1);

            uNetArray=uGrossArray-p;
            [optValue ~]=max(uNetArray,[],1);

            utilMat = reshape(optValue,nRisk, nRiskAv);

            utilMat = -exp( - RiskAvArray_2D .*  utilMat );
        end

	    % AverageCostUpdate computes the average cost for on-path contracts and computes the highest willingness to pay for off-path contracts
        function [priceTarget,isTraded] = AverageCostUpdate(Model,p)
            nRisk = length(Model.RiskSet);
            nRiskAv = length(Model.RiskAversionSet);
            nContracts = length(Model.ContractSpace);
            uGrossArray=Model.U_Gross_Array;

            uNetArray=uGrossArray-p;
            [optValue, optIndex]=max(uNetArray,[],1);            

            %Calculating average risk on-path:
            optIndexArray=repmat(optIndex,nContracts,1,1);
            contractRunVariable=repmat((1:nContracts)',1,nRisk,nRiskAv);
            
            distrArray=repmat(reshape(Model.Distribution,1,nRisk,nRiskAv),nContracts,1,1);
            riskArray=repmat(reshape(Model.RiskSet,1,nRisk,1), nContracts,1,nRiskAv);
            
            activeMass=(optIndexArray==contractRunVariable).*distrArray;
            integralRiskDistr=sum(activeMass.*riskArray,[2 3]);
            integralDistr=sum(activeMass,[2 3]);
            AC=Model.ContractSpace.*integralRiskDistr./(integralDistr+Model.BehavioralMass);
            isTraded=~isnan(AC);
            
            %Off-path contracts:
            U_OptimalArray= repmat(optValue, nContracts, 1, 1 );
            maxWillToPayArray= uGrossArray - U_OptimalArray;
            riskArray=repmat(Model.RiskSet,1,nRiskAv);
            
            riskArrayLinear=riskArray(:);
            MaxWTP_linear=reshape(maxWillToPayArray,nContracts,[]);
            
            [maxWTP, maxWTPIndex]= max(MaxWTP_linear,[], [2]);
            
            offPathPrice=min(riskArrayLinear(maxWTPIndex), maxWTP);
            
            priceTarget=AC;
            priceTarget(isnan(AC))=max(0, offPathPrice(isnan(AC)));

        end

	    % EquilibriumOperator uses AverageCostUpdate to update a guess for the price vector
        function [error,deltaPrice] = EquilibriumOperator(Model,p)
            [priceTarget,~]=AverageCostUpdate(Model,p);

            deltaPrice=priceTarget-p;
            error = norm(deltaPrice,Inf);
        end

        % DemandMatrix computes the optimal contract choice for each risk 
        % and risk aversion pair, given a price vector p
        function demandMat=DemandMatrix(model,p)
            nRisk = length(model.RiskSet);
            nRiskAv = length(model.RiskAversionSet);
            uGrossArray=model.U_Gross_Array;
            contractSpace=model.ContractSpace;

            uNetArray=uGrossArray-p;
            [~, optIndex]=max(uNetArray,[],1);

            D = reshape(optIndex,nRisk,nRiskAv);
            demandMat=contractSpace(D);
        end

        % poolMass computes the mass of consumers choosing each contract,
        % for a given price vector p
        function massInEachContract = poolMass(Model,p)
            nRisk = length(Model.RiskSet);
            nRiskAv = length(Model.RiskAversionSet);
            nContracts = length(Model.ContractSpace);
            uGrossArray=Model.U_Gross_Array;
            
            uNetArray=uGrossArray-p;
            [~, optIndex]=max(uNetArray,[],1);
            
            %Calculating mass in each contract:
            optIndexArray=repmat(optIndex,nContracts,1,1);
            contractRunVariable=repmat((1:nContracts)',1,nRisk,nRiskAv);
            distrArray=repmat(reshape(Model.Distribution,1,nRisk,nRiskAv),nContracts,1,1);
            activeMass=(optIndexArray==contractRunVariable).*distrArray;
            
            massInEachContract=sum(activeMass,[2 3]);
        end
        
        function plotTypeDistribution(modelIn, indexRiskAvLevels, save_indicator)
            arguments
                modelIn Model
                indexRiskAvLevels
                save_indicator =0
            end
%             indexRiskAvLevels = [1 5:5:30];
            if length(modelIn.RiskAversionSet) ==2
                riskAvLabels = ["\rho_l", "\rho_h"];
            else
                riskAvLabels = "\rho = " + round( modelIn.RiskAversionSet(indexRiskAvLevels), 1);
            end

            distributionNormalized = modelIn.Distribution ./ sum(modelIn.Distribution);
            lineWidth=2;
            font=15;
            fig=figure;
            hold on;
            for i=indexRiskAvLevels
                plot(modelIn.RiskSet, distributionNormalized(:,i), 'LineWidth',lineWidth);
            end
            legend( riskAvLabels, 'FontSize',font );
            title('Risk distribution conditional on \rho level', 'FontSize',font);
            xlabel('Risk level \mu', 'FontSize',font);
            ylabel('Density level', 'FontSize',font);
            set(gca,'Fontsize', 16);
            hold off;

            if save_indicator==1
                f=fullfile('../figures','DistributionPlot.pdf');
                saveas(fig,f);
            end
        end

    end
end




