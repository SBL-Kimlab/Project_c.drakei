% Gyu Min Lee: 29 nov 19
clc
dir_base = 'D:\Google_Drive_Backup\용량부족해\##Project\13.drakei_revision\Sampling_matlab\';
modelfile = 'D:\Google_Drive_Backup\용량부족해\##Project\13.drakei_revision\iSL771.mat';
genedata = 'D:\Google_Drive_Backup\용량부족해\##Project\13.drakei_revision\GIMME\Cdrakei_gene_expression.txt';
%%
% make RNA-seq data integrated model with GIMME
%integrated_model = GIMME_model(modelfile, genedata, dir_base);
%%
% do 10times of 2,000 sampling and average them
condition = ["hetero"];
array_rxn = ["wt", "FDH8", "FTHFLi", "MTHFC", "MTHFD", "MTHFR5", "METR", "CODH_ACS", "GLYCL", "GLYR"]; %wt; no change in flux(wild type)
flux_range = 0 : 0.1 : 5;

mean_mcmc_averonly(integrated_model, condition, array_rxn, flux_range, dir_base);