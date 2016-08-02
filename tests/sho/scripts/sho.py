#!/usr/bin/python

import os, sys 

if len(sys.argv)!=3:
  print "Usage: sho.py <L0> <Omega>"
  sys.exit(-1)

L0 = float(sys.argv[1])
OMEGA = float(sys.argv[2])

M=1.0/(L0**2*OMEGA)
K=OMEGA/L0**2

print M, K
