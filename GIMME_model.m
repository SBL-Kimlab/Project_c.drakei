function resultfile = GIMME_model(modelfile, genedata, dir_base)
% Make RNA-seq data integrated model.
%
% USAGE:
%   mean_sampling = mean_mcmc2(modelfile, condition, target_rxn, flux_range, dir_base)
%
% INPUTS:
%   modelfile:          COBRA model file with available format
%   genedata:           RNA-seq expression data with GeneIDs and expresiion value (TPM/FPKM/RPKM)(see mapExpressionToReactions.m)
%   dir_base:           Root of base folder you want to save the result files
%
% OUTPUTS:
%   resultfile:         full path of the integrated model
%
% EXAMPLES:
%     modelfile = 'D:\##Project\13.drakei_revision\iSL_V3.3_http_2_cobrapy.mat';
%     dir_base = 'D:\##Project\13.drakei_revision\MCMC\';
%     genedata = 'GIMME\Cdrakei_gene_expression.txt';
% 
%     integrated_model = GIMME_model(modelfile, genedata, dir_base);
%
% .. Author: - Gyu Min Lee 10/31/19

    initCobraToolbox
    changeCobraSolver('gurobi', 'all');
    [~,name,~] = fileparts(modelfile);
    if ~exist(dir_base,'dir')
        mkdir(dir_base)
        addpath(dir_base)
    end

    % load model
    model = readCbModel(modelfile);
    model = creategrRulesField(model);

    Info_genedatafile = importdata(genedata);
    expressionData.gene = Info_genedatafile.textdata;
    expressionData.value = Info_genedatafile.data(:,1);
    [expressionRxns, ~, ~] = mapExpressionToReactions(model, expressionData); 

    options.solver = 'GIMME';
    options.threshold = 0.9; %default
    options.expressionRxns = expressionRxns;

    dir_GIMME = strcat(dir_base, 'GIMME\');
    if ~exist(dir_GIMME,'dir')
        mkdir(dir_GIMME)
        addpath(dir_GIMME)
    end   
    cd(dir_GIMME)

    integrated_model = GIMME(model, options.expressionRxns, options.threshold);
    resultfile = strcat(dir_GIMME, 'GIMME_', name, '_threshold_', num2str(options.threshold), '.mat');
    if ~isfile(resultfile) 
        save(resultfile, 'integrated_model')
        fprintf ('%s ....Saved.\n\n',  resultfile);    
    else
        fprintf ('%s ....Already exist.\n\n',  resultfile);  
end    