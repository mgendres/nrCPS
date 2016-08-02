#!/usr/bin/python

import os, sys, math, time

#L = [ "8", "16", "32", "64" ]
L = [ "32" ]
#L0 = [ "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "4.0" ]
L0 = [ "3.0" ]
#OMEGA = [ "0.005", "0.010", "0.015", "0.020", "0.025", "0.030" ]
OMEGA = [ "0.005" ]
N = [ "0", "1", "2", "3", "4", "5", "6", "7" ]
#N = [ "7" ]

def get_energy(f_name):
  # Import data
  f = open(f_name, 'r')
  line = f.readline().split()
  # Prepare data for analysis
  TINY = 2.2250738585072014e-308
  dat = [ float(line[i]) + TINY for i in range(0, len(line), 2) ]
  # Compute effective mass
  effm = [ -math.log( abs( dat[i+1]/dat[i])) for i in range(0, len(dat)-1) ]
  # Determine plateau value
  count = 10
  while abs(effm[count-1]-effm[count]) - abs(effm[count+1]-effm[count]) > 0:
    if count < len(effm)-2:
#      print count, effm[count]
      count +=1
    else:
      break
  return count-1, effm[count-1]

tab = []

for l in L:
  for l0 in L0:
    for omega in OMEGA:
      for n in N:
        string_list = [ "L", str(l), "_L0", str(l0), "_Omega", str(omega), "/sho.", str(n) ]
        f_name = ''.join(string_list)
        print f_name
        count, energy = get_energy(f_name)
        theory_energy =  (float(n)+3./2.)
        if n=="7":
          theory_energy += 2.0
        theory_energy *= float(omega)
#        print "L:",l, "L0:", l0, "Omega:", omega, "N", n, "Energy:", energy/theory_energy
#        theory_energy = 1.0
        tab.append( [l, l0, omega, n, energy/theory_energy, count] )

stab = [ [ str(elem) for elem in row ] for row in tab ]

f = open("out.dat","w")
for line in stab:
  f.write(' '.join(line))
  f.write("\n")
f.close()

