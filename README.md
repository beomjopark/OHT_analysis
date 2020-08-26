# OHT_Analysis

This repository contains the reproducible code for [*Spatio-Temporal Local Interpolation For Quantifying Global Ocean Heat Transport From Autonomous Observations*](https://www.overleaf.com/read/djzmqmpsmzgn).

Data Files and resulting outputs can be accessible from NCAR GLADE file system: `./work/beomjop/OHC_dynamics/Data` and `./work/beomjop/OHC_dynamics/Results`.


## Requirements

1. The code depends on [*TEOS-10* toolbox](http://www.teos-10.org/software.htm) for `MATLAB` on your system. In the wrappers, you'll see `../gsw_matlab` which is the relative path of the TEOS-10 toolbox.


2. We assume you'll submit the job to HPC with `PBS PRO` scheduler. See `jobscript` folder for detail. 


## Pipeline

0. The `MATLAB` code are intended to call `Params_...` first to load the parameter setups (although you would need to modify some of them depending on the `PBS -v` arguments.).


### Local Interpolation

1. `selection_Int.pbs` : Profile selection and filtering procedure.

    The script consists of 4 MATLAB calls, but you should comment out all except the one you are intended to execute. Due to the computation time limit, `selectionAndVerticalIntegrationPchip_wrapper` splits the total dataset into `nParts`: Pre-2017 into 4 parts and 2017-18 into 2 parts. Check the `qsub` command to run which chunk (`iPart`) to process.

    `createDataMask_Distrib` will generate data mask for each pressure level in the list `intStartList`. `filterUsingMasks_Distrib` is actually a deprecated placeholder.


2. `fieldEM.pbs` : Mean and Covariance parameter estimation procedure.

    If you are interested in EM estimate centered at specific month(`mon`), you could just run `fieldEM.pbs` with prespecified total EM iterations(`nIterEM`). KS18 paper corresponds to `mon=2,nIterEM=0`.

     Use `fieldEM_allmonth.sh` to account for all months, i.e., `mon=1:12`, with exact EM iteration number(`iter`). This is because the all month estimation requires significantly longer time to run since we need to estimate local GP at each 12 month windows. Also, somewhat due to MATLAB license, it's often the case PBS fails time to time for certain months. Be sure to check all month estimates are computed before continuing to the next EM iteration. When running all months, be sure to match the `intStartList` for both `fieldEM_mean.pbs` and `fieldEM_cov.pbs`.

     Note that the procedures require `Distrib_Computing_Toolbox` & `Optimization_Toolbox`.


3. `anomPred_latDyn.pbs`: Kriging procedure.

    `isDeriv` specifies whether predictive derivative is of interest. When TRUE, it should accompany `targetVar` which is either 'lat' or 'lon'. Sister jobscript `anomPred_lonDyn.pbs` is intended to specify `isDeriv = true, targetVar = 'lon'`. Note that if you are doing debiasing procedure, you will see `isDeriv=false`.


### Debias

1. `adjustAnom_latDyn.pbs` : Debiasing procedure.

    This job adjusts the kriged anomaly field and prepare the divided and extended residual to refit the covariance parameters and kriging. Be sure to run `anomPred_latDyn.pbs` with `isDeriv=false` BEFORE running the debiasing procedure.


2. `fieldEM_cov_adjusted.pbs` & `anomPred_latDyn_adjusted.pbs`

    Refit covariance parameters and produce the kriged field after the correction. These codes are equivalent to scripts without `_adjusted` if `nAdjust=0` is specified.


## Misc
    
1. `integrateTKE.m`: Compute Total Kinetic Energy based on gridded velocities. This is a standalone code.



## Known Issues

1. Instability of the local GP estimate

    This requires careful investigation of (A) anomalous profiles and (B) bandwidth parameter.

2. Anomalous profiles

    Can this be ameliorated by using CORA dataset?
