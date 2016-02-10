#!/bin/bash
#PBS -N scat9
#PBS -q normal
#PBS -l nodes=8:ppn=16:native:noflash
#PBS -l walltime=00:10:00
#PBS -m abe
#PBS -M ttnhan@hotmail.com
#PBS -V
#PBS -v Catalina_maxhops=None
cd $PBS_O_WORKDIR
mkdir -p scat9
cp scat9.sh scat9/
mpirun_rsh -hostfile $PBS_NODEFILE -np 128 ./scat -n 1000000000 -p 8000 -c 27000 -k 0.5 > scat9/scat9.out 2> scat9/scat9.err
