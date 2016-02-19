#!/bin/bash
#PBS -N scat4b
#PBS -q normal
#PBS -l nodes=22:ppn=16:native:noflash
#PBS -l walltime=03:00:00
#PBS -m abe
#PBS -M ttnhan@hotmail.com
#PBS -V
#PBS -v Catalina_maxhops=None
cd $PBS_O_WORKDIR
mkdir -p scat4b
cp scat4b.sh scat4b/
ibrun --npernode 4 ./scat -n 4000000000 -p 8000 -c 27000 -k 0.5 > scat4b/scat4b.out 2> scat4b/scat4b.err
