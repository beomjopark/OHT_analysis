#!/bin/bash
# Submit the LocalMLE estimation for each month in two batches
# bash ./fieldEM_allmonth.sh
iter=3  # EM counter 0, 1, ...

## Mean
JOB_INT=$(qsub -v iter=${iter},mon=1:12 fieldEM_mean.pbs)

## Local MLE
i=1
JOB_LGP_ONE=$(qsub -W depend=afterok:${JOB_INT} -v iter=${iter},mon=${i} fieldEM_cov.pbs)
echo "${JOB_LGP_ONE}"
JOB_LGP_TWO=$(qsub -W depend=afterok:${JOB_INT} -v iter=${iter},mon=$((2*i)) fieldEM_cov.pbs)
echo "${JOB_LGP_TWO}"

JOB_INT=$JOB_LGP_ONE

for ((i=2; i<7; ++i)) ; do
    JOB_LGP_ONE=$(qsub -W depend=afterok:${JOB_INT} -v iter=${iter},mon=$((2*i-1)) fieldEM_cov.pbs)
    echo "${JOB_LGP_ONE}"

    JOB_LGP_TWO=$(qsub -W depend=afterok:${JOB_INT} -v iter=${iter},mon=$((2*i)) fieldEM_cov.pbs)
    echo "${JOB_LGP_TWO}"

    JOB_INT=$JOB_LGP_ONE
done

#qsub -v iter=${iter},mon=1:12 anomPred_latDyn.pbs
# bash ./fieldEM_allmonth.sh
