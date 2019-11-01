% Gyu Min Lee: 31 oct 19

dir_base = 'D:\Google_Drive_Backup\용량부족해\##Project\13.drakei_revision\TEST\';
modelfile = 'D:\Google_Drive_Backup\용량부족해\##Project\13.drakei_revision\iSL_V3.3_http_2_cobrapy.mat';
genedata = 'D:\Google_Drive_Backup\용량부족해\##Project\13.drakei_revision\GIMME\Cdrakei_gene_expression.txt';
%%
% make RNA-seq data integrated model with GIMME
integrated_model = GIMME_model(modelfile, genedata, dir_base);
%%
% do 10times of 2,000 sampling and average them
condition = ["auto"];
target_rxn = ["FDH7", "FTHFLi", "MTHFC", "MTHFD", "MTHFR5", "METR", "CODH_ACS", "GLYCL", "GLYR"];
flux_range = 0 : 0.1 : 5;
clc
for rxn = target_rxn
    mean_mcmc(integrated_model, condition, rxn, flux_range, dir_base);
end