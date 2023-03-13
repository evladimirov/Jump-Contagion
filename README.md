# Jump-Contagion
 
This repository illustrates and provides code for the estimation procedure developed in the paper **Jump Contagion among Stock Market Indices: Evidence from Option Markets** written by Peter Boswijk, Roger Laeven, Andrei Lalu and Evgenii Vladimirov. The paper is available on [SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3929515).

# Simulated data

The illustration is based on simulated data from the bivariate option pricing models proposed in the paper. The simulated data are stored as MATLAB .mat file in `data`. The state vectors (stored in variables mY, mV and mLambda) are simulated using an Euler discretization and the option prices (variables mOptPriceIV1 and mOptPriceIV2) are computed using the COS method. For further details concerning the simulations, see Section 3 of the paper and Appendix B.3 of Online Appendix. 

# Code

The code is written in `MATLAB` version R2019b. The `main.m` file reads the simulated options data and estimates the model parameters based on the partial-infromation C-GMM. The estimation minimizes the criterion function ‘mSVhatHJ_crit_inst4.m’, which in turn involves the implied state procedure (function ‘mSVhatHJ_ImpIntens’) and four numerical integrations of criterion functions based on the marginal states (function ‘mSVhatHJ_int_inst4’). 

The simulation and estimation procedures are described in Section 3 of the paper. 
