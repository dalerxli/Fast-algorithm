#!/bin/bash
#$ -N a_good_job_name
#$ -pe mpi-fill 16
#$ -l h_rt=00:10:00
#$ -M mail@domain.name -m bae

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
./cube 1.0 1.0 4 100 100 1.0e-6 1>${JOB_NAME}/output 2>${JOB_NAME}/errors

# move the output data to subdirectory
mv sol.dat ${JOB_NAME}
