# OHT_Analysis

This repo contains the code for "SPATIO-TEMPORAL LOCAL INTERPOLATION FOR QUANTIFYING GLOBAL OCEAN HEAT TRANSPORT FROM AUTONOMOUS OBSERVATIONS".

## Pipeline

    0.1 We assume you'll submit the job to PBS PRO scheduler. See `jobscript` folder for detail. 
    0.2 The `MATLAB` code are intended to call `Params_...` first to load the parameter setups (although you would need to modify some of them depending on the `PBS -v` arguments.).


    1. `selection_Int.pbs` : Profile selection and filtering procedure. 
    PBS consists of 4 MATLAB calls, but you should comment out all except the one you are intended to call. Due to the computation time limit, selectionAndVerticalIntegrationPchip_wrapper splits the total dataset into nParts: Pre-2017 into 4 parts and 2017-18 into 2 parts. Check the `qsub` command to run which chunk (`iPart`) to process.

    `createDataMask_Distrib` will generate data mask for each pressure level in the list `intStartList`. `filterUsingMasks_Distrib` is actually a deprecated placeholder.


    2. `fieldEM.pbs` : Mean and Covariance parameter estimation procedure.
    If you are interested in EM estimate centered at specific month(`mon`), you could just run `fieldEM.pbs` with prespecified total EM iterations(`nIterEM`). KS18 paper corresponds to `mon=2,nIterEM=1`.

    To account for all months, i.e., `mon=1:12`, use `fieldEM_allmonth.sh` with exact EM iteration number(`iter`). This is because the all month estimation requires significantly longer time to run since we need to estimate local GP at each 12 month windows. Also, somewhat due to MATLAB license, it's often the case PBS fails time to time for certain months. Be sure to check all month estimates are computed before continuing to the next EM iteration. When running allmonth, be sure to match the `intStartList` for both `fieldEM_mean.pbs` and `fieldEM_cov.pbs`.


    3. `anomPred_latDyn.pbs`: Kriging procedure.
    `isDeriv` specifies whether predictive derivative is of interest. When TRUE, it should accompany `targetVar` which is either 'lat' or 'lon'. Sister jobscript `anomPred_lonDyn.pbs` is intended to specify `isDeriv = true, targetVar = 'lon'`. Note that if you are doing debiasing procedure, you will see `isDeriv=false`.

    
    3.1 `adjustAnom_latDyn.pbs` : Debiasing procedure.


## Misc
    
    1. `integrateTKE.m`: Compute Total Kinetic Energy based on gridded velocities.
