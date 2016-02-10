#!/bin/bash
#PBS -N scat3b
#PBS -q normal
#PBS -l nodes=15:ppn=16:native:noflash
#PBS -l walltime=01:00:00
#PBS -m abe
#PBS -M ttnhan@hotmail.com
#PBS -V
#PBS -v Catalina_maxhops=None
cd $PBS_O_WORKDIR
mkdir -p scat3b
cp scat3b.sh scat3b/
mpirun_rsh -hostfile $PBS_NODEFILE -np 240 ./scat -n 3000000000 -p 8000 -c 27000 -k 0.5 -ksp_monitor -ksp_type cg -log_summary -malloc_log > scat3b/scat3b.out 2> scat3b/scat3b.err
