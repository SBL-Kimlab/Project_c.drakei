There are 3 .m files;

1. Analysis.m: In this code, you can operate GIMME_model.m and mean_mcmc.m.
2. GIMME_model.m: Using GIMME function, this code makes RNA-seq data integrated model.
3. mean_mcmc.m: There are 2 functions in this code; mcmc_10times and mean_mcmc.
   When you operate mean_mcmc function in Analysis.m, it will first run mcmc_10times
   that operate 10 times of 2,000 MCMC samplings with CHRR algorithm changing condition of model and bounds of reaction.
   After that, it will make structure with mean of MCMC sampling and name of reactions.
   A xlsx file with original sampling values and the values normalized with Biomass value. will be saved.
