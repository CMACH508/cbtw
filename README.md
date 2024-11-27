# Cauchy-Schwarz bounded trade-off weighting for causal inference with small sample sizes
The difficulty of causal inference for small-sample-size data lies in the issue of inefficiency that the variance of the estimators may be large. Some existing weighting methods adopt the idea of bias-variance trade-off, but they require manual specification of the trade-off parameters. To overcome this drawback, in this article, we propose a Cauchy-Schwarz Bounded Trade-off Weighting (CBTW) method, in which the trade-off parameter is theoretically derived to guarantee a small Mean Square Error (MSE) in estimation. We theoretically prove that optimizing the objective function of CBTW, which is the Cauchy-Schwarz upper-bound of the MSE for causal effect estimators, contributes to minimizing the MSE. Moreover, since the upper-bound consists of the variance and the squared L2-norm of covariate differences, CBTW can not only estimate the causal effects efficiently, but also keep the covariates balanced. Experimental results on both simulation data and real-world data show that the CBTW outperforms most existing methods especially under small sample size scenarios.
## Requirements
- R 4.0.5
- plyr
- ggplot2
- Rsolnp
- leaps
- MASS
- reshape2
- parallel

## Generate simulation data
```shell
Rscript sim_linear_gen.R
Rscript sim_nonlinear_gen.R
```
## Generate real-world data
```shell
Rscript ihdp_gen.R
Rscript lalonde_gen.R
```
## Run the experiments
```shell
Rscript sim_experiment.R
Rscript real_experiment.R
```
## Citation


Please cite our work if you find our code/paper is useful to your work.


```
@article{MA2025109311,
title = {Cauchy-Schwarz bounded trade-off weighting for causal inference with small sample sizes},
journal = {International Journal of Approximate Reasoning},
volume = {176},
pages = {109311},
year = {2025},
author = {Qin Ma and Shikui Tu and Lei Xu},
}
