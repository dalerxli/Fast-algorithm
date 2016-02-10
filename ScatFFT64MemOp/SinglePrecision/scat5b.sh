#!/bin/bash
#PBS -N scat5b
#PBS -q normal
#PBS -l nodes=20:ppn=16:native:noflash
#PBS -l walltime=2:00:00
#PBS -m abe
#PBS -M ttnhan@hotmail.com
#PBS -V
#PBS -v Catalina_maxhops=None
cd $PBS_O_WORKDIR
mkdir -p scat5b
cp scat5b.sh scat5b/
mpirun_rsh -hostfile $PBS_NODEFILE -np 320 ./scat -n 5000000000 -p 8000 -c 27000 -k 0.5 -ksp_monitor -log_summary -malloc_log > scat5b/scat5b.out 2> scat5b/scat5b.err
