#!/bin/bash
i="0"
while [ -f sho.$i ]
do
  ../../../../scripts/python~/corr.py sho.$i $i
i=$[$i+1]
done
