function [avgResult, xlsxFile] = averageAfterSampling(modelFile, growthCondition, rxnChangeFlux, fluxRange, defaultFolder)
% Make structure with avgerage of MCMC sampling and name of reactions and will
% make .xlsx file with original sampling values and the values normalized
% with Biomass value.
%
% USAGE:
%   [avgResult, xlsxFile] = averageAfterSampling(modelFile, growthCondition, rxnChangeFlux, fluxRange, defaultFolder)
%
% INPUTS:
%   modelFile:          COBRA model file with available format
%   growthCondition:    {'auto', 'hetero'} Condition of the model
%   rxnChangeFlux:      String array of reactions to change flux
%   fluxRange:          Array of the range of flux you want to use, if you set
%                       rxnChangeFlux as 'wt', fluxRange wll not give any influence.
%   defaultFolder:      Root of base folder you want to save the result files
%
% OUTPUTS:
%   avgResult:         Structure with mean of MCMC sampling and name of reactions
%   xlsxFile:          excel file containing name of reactions, sampling
%                      values, and sampling values normalized with biomass value
%
% EXAMPLES:
%     modelFile = 'D:\##Project\13.drakei_revision\iSL771.mat';
%     defaultFolder= 'D:\##Project\13.drakei_revision\MCMC\';
% 
%     growthCondition = ["auto"];
%     rxnChagneFlux = ["wt", "FDH8", "FTHFLi", "MTHFC", "MTHFD", "MTHFR5", "METR", "CODH_ACS", "GLYCL", "GLYR"];
%     fluxRange = 0 : 0.1 : 5;
%     clc
%     averageAfterSampling(modelFile, growthCondition, rxnChangeFlux, fluxRange, defaultFolder);
%
% .. Author: - Gyu Min Lee 06/17/20

    initCobraToolbox;
    changeCobraSolver('ibm_cplex', 'all');
    [~,modelFileName,~] = fileparts(modelFile);
    
    if ~exist(defaultFolder,'dir')
        mkdir(defaultFolder)
        addpath(defaultFolder)
    end

    for rxn = rxnChangeFlux
        if strcmp(growthCondition, 'auto') | strcmp(growthCondition, 'hetero')   %#ok<OR2>
            growthCondition1 = strcat(growthCondition, 'trophic');
        else
            growthCondition1 = '';
        end
        fprintf ('%s\n', growthCondition1);
        growthConditionFolder = strcat(defaultFolder, growthCondition1);
        if ~exist(growthConditionFolder, 'dir')
            mkdir(growthConditionFolder)
            addpath(growthConditionFolder)
        end
        resultFolder = strcat(growthConditionFolder, '\', rxn);        
        if ~exist(resultFolder, 'dir')
            mkdir(resultFolder)
            addpath(resultFolder)
        end
        cd(resultFolder)
           
        for flux = fluxRange
            if ~strcmp(rxn, 'wt')
                fprintf ('%s, %s flux from %s to %s ....operating\n\n', modelFileName, rxn, num2str(fluxRange(1,1)), num2str(fluxRange(1,end)));
            elseif strcmp(rxn, 'wt')
                fprintf ('%s, %s ....operating\n\n', modelFileName, rxn);
            end
            [modelChange, samplingResult] = mcmc(modelFile, growthCondition, rxn, flux);
            
            [~,samplingResultFileName,~] = fileparts(samplingResult);
            model = strsplit(samplingResultFileName, '_');
        
            if ~isfile(samplingResult)
                if ~strcmp(rxn, 'wt')
                    warning (sprintf ('%s : flux %s .... NOT finished.\nend',  rxn, num2str(flux)));
                elseif strcmp(rxn, 'wt')
                    warning (sprintf ('%s .... NOT finished.\nend',  rxn));
                end
            else
                if ~strcmp(rxn, 'wt')
                    fprintf('%s : flux %s .... finished\nNow making average ....\n',  rxn, num2str(flux))
                elseif strcmp(rxn, 'wt')
                    fprintf('%s .... finished\nNow making average ....\n',  rxn)
                end
                avgFile = strcat('avg_', samplingResultFileName, '_100000.mat');
                xlsxFile = strcat('avg_', samplingResultFileName, '_100000.xlsx');
                if ~isfile (avgFile)
                    load (samplingResult);
                    avgResult = mean(samples, 2);
                    avgResult = struct('avg', avgResult);
                    [avgResult(:).rxns] = modelChange.rxns;
                    save(avgFile, 'avgResult');
                    fprintf ('%s ....Saved.\n',  avgFile);

                    avgResultSize = size(avgResult.avg, 1);
                    avgNormalize = zeros(avgResultSize, 1);
                    idxBiomass = find(modelChange.c);
                    valueBiomass = avgResult.avg(idxBiomass);
                    for i = 1 : avgResultSize
                        avgNormalize(i, 1) = avgResult.avg(i, 1) / valueBiomass;
                        avgNormalize(i, 1) = round(avgNormalize(i, 1), 5);
                    end
                    mergeRxnAndAvg = [avgResult.rxns, num2cell(avgResult.avg), num2cell(avgNormalize)];
                    writecell(mergeRxnAndAvg, xlsxFile, 'Sheet', 1, 'Range', 'A1')
                    fprintf ('%s ....Saved.\n\n',  xlsxFile);
                else
                    fprintf ('%s ....Already exist.\n',  avgFile);
                    if ~isfile(xlsxFile)
                        load(avgFile);

                    avgResultSize = size(avgResult.avg, 1);
                    avgNormalize = zeros(avgResultSize, 1);
                    idxBiomass = find(modelChange.c);
                    valueBiomass = avgResult.avg(idxBiomass);
                    for i = 1 : avgResultSize
                        avgNormalize(i, 1) = avgResult.avg(i, 1) / valueBiomass;
                        avgNormalize(i, 1) = round(avgNormalize(i, 1), 5);
                    end
                    mergeRxnAndAvg = [avgResult.rxns, num2cell(avgResult.avg), num2cell(avgNormalize)];
                    writecell(mergeRxnAndAvg, xlsxFile, 'Sheet', 1, 'Range', 'A1')
                    fprintf ('%s ....Saved.\n\n',  xlsxFile);                 
                    else
                        fprintf ('%s ....Already exist.\n\n',  xlsxFile);
                    end           
                end    
            end
        end
    end
end

function [modelChange, samplingResult] = mcmc(modelFile, growthCondition, rxn, flux)
% Operate 100,000 MCMC sampling with CHRR algorithm changing
% condition of model and bounds of reaction.
%
% USAGE:
%   [modelChange, samplingResult] = mcmc(modelFile, growthCondition, rxn, flux)
%
% INPUTS:
%   modelFile:          COBRA model file with available format
%   growthCondition:    {'auto', 'hetero'} Condition of the model
%   rxn:                String of reaction to change flux (wild type; 'wt')
%   flux:               Double of flux you want to use
%   defaultFolder:      Root of base folder you want to save the result files
%
% OUTPUTS:
%   modelChange:        Flux changed model used in smapling, providing a lower bound of
%                       95% of the optimal growth rate as computed by FBA
%   samplingResult:     Double with 100,000 times of MCMC sampling
%
% .. Author: - Gyu Min Lee 06/17/20
    
    [~,modelFileName,~] = fileparts(modelFile);
    model = readCbModel(modelFile);
    objectiveCol = [find(model.c)];
    options.nPointsReturned = 100000;

    switch growthCondition
        case 'autotrophic'
            model = changeRxnBounds( model, 'EX_h2_e', -10, 'l' );
            model = changeRxnBounds( model, 'EX_co2_e', -5, 'l' );
            model = changeRxnBounds( model, 'EX_fru_e', 0, 'l' );
        case 'heterotrophic'
            model = changeRxnBounds( model, 'EX_h2_e', 0, 'l' );
            model = changeRxnBounds( model, 'EX_co2_e', 0, 'l' );
            model = changeRxnBounds( model, 'EX_fru_e', -2.224, 'l' );
        otherwise
            model = model;         
    end
    
    modelChange = model;
    if ~strcmp(rxn, 'wt')
        samplingResult = strcat(modelFileName, '_', rxn, '_flux_', num2str(flux), '_100000sampling', '.mat');
        modelChange = changeRxnBounds( modelChange, rxn, flux, 'b');
    elseif strcmp(rxn, 'wt')
        samplingResult = strcat(modelFileName, '_', rxn, '_100000sampling', '.mat'); 
    end
    
    if ~isfile(samplingResult)    
        fba = optimizeCbModel(modelChange);
        modelChange = changeRxnBounds(modelChange, modelChange.rxns(objectiveCol, 1), 0.95 * fba.f, 'l');                           
        try
            warning('off')
            [ ~,samples ] = sampleCbModel( modelChange,[], [], options );
            warning('on')
            if ~strcmp(rxn, 'wt')
                fprintf ('%s : flux %s, %s ....operate.\n',  rxn, num2str(flux), '100,000 times');
            elseif strcmp(rxn, 'wt')
                fprintf ('%s : %s ....operate.\n',  rxn, '100,000 times');
            end
            save(samplingResult, 'samples')
            fprintf ('%s ....Saved.\n\n',  samplingResult);    
        catch
            warning('on')
            if ~strcmp(rxn, 'wt')
                warning (sprintf ('%s : flux %s, %s ....NOT operate.',  rxn, num2str(flux), '100,000 times'));
            elseif strcmp(rxn, 'wt')
                warning (sprintf ('%s : %s ....NOT operate.',  rxn, '100,000 times'));
            end
        end
    else
        fprintf ('%s ....Already exist.\n\n',  samplingResult);            
    end
end    