#!/bin/bash
#PBS -N scat4b
#PBS -q normal
#PBS -l nodes=16:ppn=16:native:noflash
#PBS -l walltime=01:00:00
#PBS -m abe
#PBS -M ttnhan@hotmail.com
#PBS -V
#PBS -v Catalina_maxhops=None
cd $PBS_O_WORKDIR
mkdir -p scat4b
cp scat4b.sh scat4b/
mpirun_rsh -hostfile $PBS_NODEFILE -np 256 ./scat -n 4000000000 -p 8000 -c 27000 -k 0.5 -ksp_monitor -ksp_type cg -log_summary -malloc_log > scat4b/scat4b.out 2> scat4b/scat4b.err
