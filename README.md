There are 3 .m files;

1. Analysis.m;
   In this code, you can operate OperateGIMME.m and averageAfterSampling.m.

2. OperateGIMME.m;
   Using GIMME function, this code makes RNA-seq data integrated model.

3. averageAfterSampling.m:
   There are 2 functions in this code; mcmc and averageAfterSampling. 
   When you operate averageAfterSampling function in Analysis.m, it will first run 
   mcmc that operate 100,000 MCMC samplings with CHRR algorithm changing condition
   of model and bounds of reaction.   
   After that, it will make structure with average of MCMC sampling and name of reactions.
   A xlsx file with original sampling values and the values normalized with Biomass value will be saved.
