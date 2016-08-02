#!/bin/bash
i="0"
while [ -f slater.$i.$i ]
do
  ../../../../scripts/python/corr2.py slater.$i.$i $i
#  ../../../../scripts/python/signoise.py slater.$i $i
i=$[$i+1]
done
