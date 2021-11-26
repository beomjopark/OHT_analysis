#!/bin/bash
### Job Name
#PBS -N Anom_Pred_Dyn_DUACS
### Project code
#PBS -A UMEL0001
#PBS -l walltime=12:00:00
#PBS -q economy
### Merge output and error files
#PBS -j oe
### Write stdout and stderr on runtime
#PBS -k oe
### Select 1 node with 36 CPUs
#PBS -l select=1:ncpus=36
### Send email on abort, begin and end
#PBS -m abe
### Specify mail recipient
#PBS -M beomjop@andrew.cmu.edu

##$ Load Matlab
module load matlab/R2020a

### Use bash
export SHELL=/bin/bash

### Run the executable
echo "iterEM: ${iter}"
echo "Month: ${mon}"
echo "Type: ${type}"
echo "isProfile: ${isProfile}"
echo "isAdjusted: ${isAdjusted}"

matCMD="Params_LatFlux_Step1, isAdjusted = ${isAdjusted}, "
if [ ${isAdjusted} = true ]
then
    echo "Adjusted"
    matCMD+="nAdjust=1, "
fi

matCMD+="isDeriv = true, targetVar = '${type}', "
matCMD+="month = ${mon}, "
matCMD+="iterEM = ${iter}, "
matCMD+="isProfile = ${isProfile}, "
matCMD+="isResetRes = false, "
matCMD+="responseTag = 'DUACS', "

if [ ${mon} = 1:12 ]
then
    echo "FullMonth!"
    matCMD+="isFullMonth = true, "
fi
#intStartList = [10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800] % Full upper ocean
matCMD+="intStartList = [ 10], "
matCMD+="computeAnomalies_wrapper, "
matCMD+="exit"
echo "${matCMD}"

matlab -nodisplay -r "${matCMD}"
# Lat Flux
#matlab -nodisplay -r "Params_LatFlux_Step1, isAdjusted = false, computeAnomalies_wrapper, exit"

# qsub -v iter=3,mon=11 anomPred_latDyn.pbs