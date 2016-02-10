#!/bin/bash
#$ -N Scattering3DS-15
#$ -pe mpi-fill 80
#$ -l h_rt=24:10:00
#$ -M nhantran@ksu.edu -m bae

###############################################################################
# Notes for Sun Grid Engine Lines
###############################################################################
# 
# 1) Change the line '-N ...' to a reasonable job name.
# 2) Change the line '-pe ...' to use a good number of CPUs.
# 3) Change the line '-l h_rt ..." to a good overestimate of the runtime.
# 4) Change the line '-M ...' to use the correct email address.
#
###############################################################################


# change directory to working directory
cd ${SGE_O_WORKDIR}

# make a subdirectory to place data and output in
mkdir -p ${JOB_NAME}

# copy the job script to the subdirectory
cp run.sh ${JOB_NAME}

# run the simulation (NOTE: change parameters here)
./scattering3DS 1.0e-15 10000 10000 1.0e-6 1>${JOB_NAME}/output 2>${JOB_NAME}/errors

# move the output data to subdirectory
mv sol.dat ${JOB_NAME}
