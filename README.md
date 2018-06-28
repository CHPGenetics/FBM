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
