% Gyu Min Lee: 17 june 20
clc
defaultFolder = 'D:\Google_Drive_Backup\용량부족해\##Project\13.drakei_revision\Sampling_matlab\t22\';
modelFile = 'D:\Google_Drive_Backup\용량부족해\##Article\project_drakei\3. 2019_PNAS_yoseb\1.before_revision\input\iSL771.mat';
geneData = 'D:\Google_Drive_Backup\용량부족해\##Project\13.drakei_revision\GIMME\Cdrakei_gene_expression.txt';
%%
% make RNA-seq data integrated model with GIMME
integratedModelFile = OperateGIMME(modelFile, geneData, defaltFolder);
%%
% do 100,000 sampling and average them
growthCondition = ["auto", "hetero"]; % auto, hetero, or empty string("")
rxnChangeFlux = ["wt", "FDH8", "FTHFLi", "MTHFC", "MTHFD", "MTHFR5", "METR", "CODH_ACS", "GLYCL", "GLYR"]; %wt; no change in flux(wild type)
fluxRange = 0 : 0.1 : 5;

averageAfterSampling(integratedModelFile, growthCondition, rxnChangeFlux, fluxRange, defaultFolder);