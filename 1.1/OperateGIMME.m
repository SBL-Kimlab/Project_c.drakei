function integratedModelFile = OperateGIMME(modelFile, geneData, defaltFolder)
% Make RNA-seq data integrated model.
%
% USAGE:
%   integratedModelFile = OperateGIMME(modelFile, geneData, defaltFolder)
%
% INPUTS:
%   modelFile:            COBRA model file with available format
%   geneData:             RNA-seq expression data with GeneIDs and expresiion value (TPM/FPKM/RPKM)(see mapExpressionToReactions.m)
%   defaultFolder:        Root of base folder you want to save the result files
%
% OUTPUTS:
%   integratedModelFile:  full path of the integrated model
%
% EXAMPLES:
%     modelFile = 'D:\##Project\13.drakei_revision\iSL771.mat';
%     defaultFolder = 'D:\##Project\13.drakei_revision\MCMC\';
%     geneData = 'GIMME\Cdrakei_gene_expression.txt';
% 
%     integratedModel = OperateGIMME(modelFile, geneData, defaultFolder);
%
% .. Author: - Gyu Min Lee 06/17/20

    initCobraToolbox
    changeCobraSolver('gurobi', 'all');
    [~,name,~] = fileparts(modelFile);
    
    if ~exist(defaltFolder,'dir')
        mkdir(defaltFolder)
        addpath(defaltFolder)
    end

    % load model
    model = readCbModel(modelFile);
    model = creategrRulesField(model);

    geneDataInfo = importdata(geneData);
    expressionData.gene = geneDataInfo.textdata;
    expressionData.value = geneDataInfo.data(:,1);
    [expressionRxns, ~, ~] = mapExpressionToReactions(model, expressionData); 

    options.solver = 'GIMME';
    options.threshold = 0.9; %default
    options.expressionRxns = expressionRxns;

    GIMMEFolder = strcat(defaltFolder, 'GIMME\');
    if ~exist(GIMMEFolder,'dir')
        mkdir(GIMMEFolder)
        addpath(GIMMEFolder)
    end   
    cd(GIMMEFolder)

    integratedModel = GIMME(model, options.expressionRxns, options.threshold);
    integratedModelFile = strcat(GIMMEFolder, 'GIMME_', name, '_threshold_', num2str(options.threshold), '.mat');
    if ~isfile(integratedModelFile) 
        save(integratedModelFile, 'integratedModel')
        fprintf ('%s ....Saved.\n\n',  integratedModelFile);    
    else
        fprintf ('%s ....Already exist.\n\n',  integratedModelFile);
    end
    cd (defaltFolder)
end    