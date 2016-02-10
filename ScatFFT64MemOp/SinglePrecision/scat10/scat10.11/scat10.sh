#!/bin/bash
#PBS -N scat10
#PBS -q normal
#PBS -l nodes=40:ppn=16:native:noflash
#PBS -l walltime=5:00:00
#PBS -m abe
#PBS -M ttnhan@hotmail.com
#PBS -V
#PBS -v Catalina_maxhops=None
cd $PBS_O_WORKDIR
mkdir -p scat10
cp scat10.sh scat10/
ibrun -N 40 -npernode 5 ./scat -n 10000000000 -p 8000 -c 27000 -k 0.5 -ksp_monitor -log_summary -malloc_log > scat10/scat10.out 2> scat10/scat10.err
