# BCVAR
Code that replicates the Bayesian Compressed Vector Autoregressive (BCVAR) model in Koop, G., Korobilis, D. and Pettenuzzo, D. (2019). “Bayesian Compressed Vector Autoregressions”, Journal of Econometrics, 210, 135-154. 

This code is for demonstration of the Bayesian compressed VAR (BCVAR) method. Because it is non-trivial (computationally) to replicate all our results and extensive robustness checks, we include a main function FORECASTING.m set up to forecast using the MEDIUM-sized VARs for the various methods we use in our paper. Forecasts also depend on the posterior mean of the predictive density (no predictive simulation due to requirements in RAM), that is, there is no simulation. The default settings should be accessible to a wide-range of researchers interested in understanding how BCVAR works.

If you are only interested in understanding how compressed VAR works, then you can examine the MATLAB functions
…\functions\BCTRVAR_CONJ.m     % estimates triangular BCVAR, with option to compress cov matrix
…\functions\BCTRVAR_TVP.m        % estimates TVP version of the triangular BCVAR
We also have a function that estimates BCVAR using the natural conjugate prior. As we explain in the paper, this case does not work that well and it is very restrictive since it restricts VAR equations to be symmetric (same variables showing up in all VAR equations). This is the MATLAB function
…\functions\BCVAR_CONJ.m          % estimates BCVAR using natural conjugate prior

It is straightforward to ask the code to estimate larger models and also do predictive simulation (so that the full predictive density is available, not just its mean). For example, in the case of the LARGE VAR all is needed is to set:
VAR_size = 'LARGE';
ndraws = 10000;
Nevertheless, these settings will require a large amount of resources (CPU and RAM).  What we did in practice was to split estimation over different methods/models and forecast horizons in different cores of a High Performance Cluster. Therefore, attempting more demanding models is only advised for more experienced users. We provide no support on how to do this, the code is only for demonstration and it is the user’s responsibility to understand whether it suits their needs. While the default code we provide is very simple (in order to be accessible to less experienced users), more advanced programmers will be able to find all the functions we used to produce our results in the folder “functions”. Feel free to modify these functions as you wish.
