#!/bin/bash
#PBS -N scat6b
#PBS -q normal
#PBS -l nodes=40:ppn=16:native:noflash
#PBS -l walltime=02:00:00
#PBS -m abe
#PBS -M ttnhan@hotmail.com
#PBS -V
#PBS -v Catalina_maxhops=None
cd $PBS_O_WORKDIR
mkdir -p scat6b
cp scat6b.sh scat6b/
ibrun --npernode 8 ./scat -n 6000000000 -p 8000 -c 27000 -k 0.5 > scat6b/scat6b.out 2> scat6b/scat6b.err
