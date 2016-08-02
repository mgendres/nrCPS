#!/usr/bin/python

import os, sys 

if len(sys.argv)!=2:
  print "Usage: profile.py <file>"
  sys.exit(-1)
FILE = sys.argv[1]

if not os.path.exists(FILE):
  print ""
  print "File does not exist, exiting program."
  sys.exit(-1)

print "Parsing file...",
file = open(FILE,"r")

function_types = ["main","Lattice::Refresh","Propagator::Run","TwoBody::Run","SlaterDet2::Run"]
function_type = []
function_cost = []

for line in file:
  split_line = line.split()
  try:
    for function in function_types:
      if split_line[1][4:len(function)+4]==function:
        cost=int(split_line[0].replace(",",""))
#        print function, cost
        function_type.append( function )
        function_cost.append( int(cost) )
  except IndexError:
    pass
file.close()
print "done."

main_cost=function_cost[0]
simulation_cost=sum(function_cost[1:])
print ""
for i in range(len(function_cost)):
  print "Total time spent on", function_type[i], ":", function_cost[i]

print ""
print "Percent time spent on initialization overhead:", 1.0-float(simulation_cost)/main_cost

print ""
for i in range(1,len(function_cost)):
  print "Percent time spent on", function_type[i], ":", float(function_cost[i])/main_cost

print ""
for i in range(1,len(function_cost)):
  print "Percent time excluding overhead spent on", function_type[i], ":", float(function_cost[i])/simulation_cost

