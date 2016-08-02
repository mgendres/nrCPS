#!/bin/bash
#
# A script to perform omega and L0 scans for sho test program.
# This script should be run within the simulation directory.
# For each omega value, L0, M and K will be updated appropriately
# and output from the simulation saved in a directory specified by L0 and omega

CURRENTPATH=$PWD

L=" 8 16 32 64"
L0=" 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0"
OMEGA="0.005 0.010 0.015 0.020 0.025 0.030"

L="32"
L0=" 3.0"
OMEGA="0.005"

for l in $L
do

for l0 in $L0
do

  for omega in $OMEGA
  do
    # Determine M and K...
    echo L0: $l0 OMEGA: $omega
    PARAMS=(`$CURRENTPATH/scripts/sho.py $l0 $omega`)
    M=${PARAMS[0]}
    K=${PARAMS[1]}
    echo M: $M K: $K

    # Update appropriate arg files..
    sed -i "2c\\$M" $PWD/args/kinetic.arg
    sed -i "3c\\$K" $PWD/args/potential.arg
    sed -i "4c\\$K" $PWD/args/potential.arg
    sed -i "5c\\$K" $PWD/args/potential.arg
    sed -i "2c\\$l0" $PWD/args/one_body.arg
    sed -i "3c\\$l0" $PWD/args/one_body.arg
    sed -i "4c\\$l0" $PWD/args/one_body.arg
    sed -i "2c\\$l" "$PWD/args/do.arg"
    sed -i "3c\\$l" "$PWD/args/do.arg"
    sed -i "4c\\$l" "$PWD/args/do.arg"



    # Run...
    mkdir -p $CURRENTPATH/results
    rm $CURRENTPATH/results/*
    $CURRENTPATH/a.out

    # Compute effective masses
    cd $CURRENTPATH/results
    $CURRENTPATH/scripts/run.sh
    cd $CURRENTPATH

    # Move results directory...
    RESULTS=L${l}_L0${l0}_Omega${omega}
    mv $CURRENTPATH/results $CURRENTPATH/$RESULTS


  done
done
done

mkdir -p $PWD/results


