#!/bin/bash
#PBS -N scat10a
#PBS -q normal
#PBS -l nodes=35:ppn=16:native:noflash
#PBS -l walltime=02:00:00
#PBS -m abe
#PBS -M ttnhan@hotmail.com
#PBS -V
#PBS -v Catalina_maxhops=None
cd $PBS_O_WORKDIR
mkdir -p scat10a
cp scat10a.sh scat10a/
mpirun_rsh -hostfile $PBS_NODEFILE -np 560 ./scat -n 10000000000 -p 8000 -c 27000 -k 0.5 > scat10a/scat10a.out 2> scat10a/scat10a.err
