#  averageAfterSampling

Matlab code files to get average of 100,000 MCMC samplings with CHRR algorithm

## Description of Matlab code files

`1. Analysis.m`

   --> In this code, you can operate `2. OperateGIMME.m` and `3. averageAfterSampling.m`.

`2. OperateGIMME.m`

   --> Using GIMME function, this code makes RNA-seq data integrated model.

`3. averageAfterSampling.m`

   --> There are 2 functions in this code; mcmc and averageAfterSampling.
   When you operate averageAfterSampling function in `1. Analysis.m`, it will first run 
   mcmc that operate 100,000 MCMC samplings with CHRR algorithm changing condition 
   of model and bounds of reaction.
   After that, it will make structure taking average of MCMC sampling and name of reactions.
   A xlsx file with original sampling values and the values normalized with Biomass value will be saved.
   

## Citation

Song, Y., Lee, J. S., Shin, J., Lee, G. M., Jin, S., Kang, S., ... & Cho, S. (2020). Functional cooperation of the glycine synthase-reductase and Wood–Ljungdahl pathways for autotrophic growth of Clostridium drakei. Proceedings of the National Academy of Sciences, 117(13), 7516-7523.
