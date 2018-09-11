# FBM

## Basic user manual 
In R, call *BI.cpp* using 
```
> sourceCpp(paste(PATH, "/FBM.cpp", sep=""))
```
After data is stored as a list `data`, following instructions in *example.R*, call full MCMC run by
```
> results <- FBMcpp(data, nItr, seed, "BI")
```
MCMC results are saved in list `results`, to be estimated. Despite what's provided in this repo, different versions of estimating and plotting code can be applied to `results` tailored to users' specific need.

*In case issues raised without being fixed in a timely matter, please contact author directly at fangz.ark@gmail.com*

_Reference_ Zhou Fang, Tianzhou Ma, Gong Tang, Li Zhu, Qi Yan, Ting Wang, Juan C Celedón, Wei Chen, George C Tseng; Bayesian integrative model for multi-omics data with missingness, Bioinformatics, , bty775, https://doi.org/10.1093/bioinformatics/bty775
