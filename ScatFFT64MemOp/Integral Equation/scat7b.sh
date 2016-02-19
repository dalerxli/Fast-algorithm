#!/bin/bash
#PBS -N scat7b
#PBS -q normal
#PBS -l nodes=55:ppn=16:native:noflash
#PBS -l walltime=02:00:00
#PBS -m abe
#PBS -M ttnhan@hotmail.com
#PBS -V
#PBS -v Catalina_maxhops=None
cd $PBS_O_WORKDIR
mkdir -p scat7b
cp scat7b.sh scat7b/
ibrun --npernode 5 ./scat -n 7000000000 -p 8000 -c 64000 -k 0.5 -a 0.00001 -ksp_rtol 0.001 > scat7b/scat7b.out 2> scat7b/scat7b.err
