# SafeLabor
Joint Modeling of Competing Risks and Current Status Data: An Application to Spontaneous Labor Study.

## Code
- [ReadData.R](https://github.com/youjin1207/SafeLabor/blob/master/Code/ReadData.R) : Read the Consortium on Safe Labor Data.
- [Lognormal_sim.R](https://github.com/youjin1207/SafeLabor/blob/master/Code/Lognormal_sim.R) : One example of generating simulation data with log-normal frailty model and estimating parameters using maximum likelihood method. 
- [Weibull_sim.R](https://github.com/youjin1207/SafeLabor/blob/master/Code/Weibull_sim.R) : One example of generating simulation data with Weibull frailty model and estimating parameters using maximum likelihood method. 
- Auxiliary functions : [real_splines.R](https://github.com/youjin1207/SafeLabor/blob/master/Code/real_splines.R) \& [sim_splines.R](https://github.com/youjin1207/SafeLabor/blob/master/Code/sim_splines.R). These two are slightly different because of different starting point.
- [sim_model.R](https://github.com/youjin1207/SafeLabor/blob/master/Code/sim_model.R) : generates the explanatory plots for simulation scenarios. 
- [Lognormal_real.R](https://github.com/youjin1207/SafeLabor/blob/master/Code/Lognormal_real.R) : Based on the real data from nulliparous women and multiparous women, fit the joint model and estimate the parameters including the effect of baseline covariates. We assume three different Cox model depending on baseline covariates and assume log-normal distribution for frailty model.
- [Weibull_real.R](https://github.com/youjin1207/SafeLabor/blob/master/Code/Weibull_real.R) : Based on the real data from nulliparous women and multiparous women, fit the joint model and estimate the parameters including the effect of baseline covariates. We assume three different Cox model depending on baseline covariates and assume Weibull distribution for frailty model.
- [gradient.R](https://github.com/youjin1207/SafeLabor/blob/master/Code/gradient.R) : Diagnostic gradient function for modal diagnostics.
- Read result : [Lognormal_readresult.R](https://github.com/youjin1207/SafeLabor/blob/master/Code/Lognormal_readresult.R) \& [Weibull_readresult.R](https://github.com/youjin1207/SafeLabor/blob/master/Code/Weibull_readresult.R)
- [Rcpp folder](https://github.com/youjin1207/SafeLabor/tree/master/Code/Rcpp) : contains c++ file for calculating likelihood function used for `optim` in `R`. Using c++ (Rcpp) substantially boosts the speed. 

## Data
- [Seed folder](https://github.com/youjin1207/SafeLabor/tree/master/Data/Seed) : includes the seed information used for generating bootstrap samples for rel data.
- [lognormal_real folder](https://github.com/youjin1207/SafeLabor/tree/master/Data/lognormal_real) : contains (bootstrap) estimates from each join model assuming log-normal distribution for frailty model.
- [Weibull_real folder](https://github.com/youjin1207/SafeLabor/tree/master/Data/Weibull_real) : contains (bootstrap) estimates from each join model assuming Weibull distribution for frailty model.

## Figure
- Explanatory figures for real data : [NullWithout](https://github.com/youjin1207/SafeLabor/blob/master/Figure/NullWithout.pdf), [MultiWithout](https://github.com/youjin1207/SafeLabor/blob/master/Figure/MultiWithout.pdf), [null_hist](https://github.com/youjin1207/SafeLabor/blob/master/Figure/null_hist.pdf), \& [multi_hist_beyond10](https://github.com/youjin1207/SafeLabor/blob/master/Figure/multi_hist_beyond10.pdf).
- Explanatory figures for simulation study : [baselines](https://github.com/youjin1207/SafeLabor/blob/master/Figure/baselines.pdf) \& [frailty](https://github.com/youjin1207/SafeLabor/blob/master/Figure/frailty.pdf).
- Real data results
