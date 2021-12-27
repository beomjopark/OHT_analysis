#!/bin/bash
### Job Name
#PBS -N EMFields_intDens_NoTrend
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

matCMD="Params_LatFlux_Step1, "
matCMD+="month = ${mon}, "
matCMD+="nIterEM = ${iter}, "
matCMD+="responseTag = 'DUACS', "

if [ ${mon} = 1:12 ]
then
    echo "FullMonth!"
    matCMD+="isFullMonth = true, "
fi

matCMD+="intStartList = [ 10], "
matCMD+="estimateFieldEM_wrapper, exit"
echo "${matCMD}"

matlab -nodisplay -r "${matCMD}"

# qsub -v mon=11,iter=3 fieldEM.pbs
