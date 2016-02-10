#!/bin/bash
#PBS -N scat2b
#PBS -q normal
#PBS -l nodes=15:ppn=16:native:noflash
#PBS -l walltime=01:00:00
#PBS -m abe
#PBS -M ttnhan@hotmail.com
#PBS -V
#PBS -v Catalina_maxhops=None
cd $PBS_O_WORKDIR
mkdir -p scat2b
cp scat2b.sh scat2b/
mpirun_rsh -hostfile $PBS_NODEFILE -np 240 ./scat -n 2000000000 -p 8000 -c 27000 -k 0.5 -ksp_monitor -ksp_type cg -log_summary -malloc_log > scat2b/scat2b.out 2> scat2b/scat2b.err
