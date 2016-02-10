#!/bin/bash
#$ -N scat9c
#$ -pe mpi-fill 16
#$ -l h_rt=24:00:00
#$ -M nhantran@ksu.edu -m bae
#$ -l mem=45G
#$ -cwd
# -o output


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
cp ${JOB_NAME}.sh ${JOB_NAME}

# run the simulation (NOTE: change parameters here)
mpirun -np $NSLOTS ./scat -n 1000000000 -k 0.5 -p 8000 -c 27000 -ksp_type fgmres -ksp_monitor -log_summary >${JOB_NAME}/output1 2>${JOB_NAME}/error1

# move the output data to subdirectory
# mv sol.dat ${JOB_NAME}
