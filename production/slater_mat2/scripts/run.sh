#!/bin/bash

N=7

for i in $(seq 0 $N)
do

  for j in $(seq 0 $N)
  do

    ../../../../scripts/python/corr.py slater.$i.$j $i.$j
  done
done
