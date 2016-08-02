#!/bin/csh
qsub -I -V -q dirac_int -l nodes=1:ppn=1:tesla
