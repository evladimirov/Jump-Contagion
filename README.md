# Jump-Contagion
 
This repository illustrates and provides code for the estimation procedure developed in the paper **Jump Contagion among Stock Market Indices: Evidence from Option Markets** written by Peter Boswijk, Roger Laeven, Andrei Lalu and Evgenii Vladimirov. The latest version of the paper is available via this [link](https://www.evladimirov.com/files/BLLV_JumpContagion.pdf).

# Simulated data

The illustration is based on simulated data from the bivariate option pricing models proposed in the paper. The simulated data are stored as the MATLAB .mat file `sim.mat` in the folder `data`. The state vectors (stored in the variables mY, mV and mLambda) are simulated using an Euler discretization and the option prices (variables mOptPriceIV1 and mOptPriceIV2) are computed using the COS method. For further details concerning the simulations, see Section 4 of the paper and Appendix C.4 of the Online Appendix. 

# Code

The code is written in `MATLAB` version R2018b. The main building blocks of the `main.m` file read the simulated options data and estimate the parameters of the bivariate model proposed in the paper based on the partial-information C-GMM procedure. The estimation procedure minimizes the criterion function ‘./code/mSVhatHJ_crit_inst4.m’, which in turn involves the implied state procedure (function ‘./code/mSVhatHJ_ImpIntens.m’) and four numerical integrations of criterion functions based on the marginal states (function ‘./code/mSVhatHJ_int_inst4.m’). Given the estimated parameters, the standard errors are calculated using the function ‘./code/mSVhatHJ_std4.m’. The estimated parameters and the figure with implied intensities are displayed as the result of the optimization. The simulation and estimation procedures are described in Section 4 of the paper. 

# Output

We summarize all the replication results, including Table 3 (bivariate parameter estimates), Figure 2 (implied intensities), Table 4 (option pricing fit), Table 5 (univariate parameter estimates), Table 6 (descriptive statistics log-return distribution), Figure 3 (contour plots), Table 2 (Monte Carlo simulation results - summary) and Table C.1 (Monte Carlo simulation results - full details) from the paper, and the function dependences in Section 4 of the accompanying [pdf file](https://github.com/evladimirov/Jump-Contagion/blob/main/replication_results.pdf). 

