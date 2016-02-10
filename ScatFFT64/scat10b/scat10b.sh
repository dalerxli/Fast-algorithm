#!/bin/bash
#PBS -N scat10b
#PBS -q normal
#PBS -l nodes=64:ppn=16:native:noflash
#PBS -l walltime=10:00:00
#PBS -m abe
#PBS -M ttnhan@hotmail.com
#PBS -V
#PBS -v Catalina_maxhops=None
cd $PBS_O_WORKDIR
mpirun_rsh -hostfile $PBS_NODEFILE -np 1024 ./scat -n 10000000000 -p 8000 -c 27000 -k 0.5 -log_summary
