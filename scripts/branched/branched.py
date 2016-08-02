#!/usr/bin/python

import os, sys 

def cost(levels, steps, branches):
  return steps*branches*(branches**levels - 1)/(branches - 1)

# Branch parameters assuming T = levels*steps = 96
# level, steps, branches (per step)
PARAMS=[

[2,48,2],
[2,48,3],
[2,48,4],
[2,48,5],
[2,48,6],

[4,24,2],
[4,24,3],
[4,24,4],
[4,24,5],
[4,24,6],

[6,16,2],
[6,16,3],
[6,16,4],
[6,16,5],

[8,12,2],
[8,12,3],

[12,8,2],

]

# Compute the cost of each parameter set
COSTS=[]
for params in PARAMS:
  COSTS.append( cost(params[0], params[1], params[2]) )

# Determine the most costly parameter set
max_cost = max(COSTS)

# Compute number of configurations per block for each parameter set, assuming a fixed cost (to within 50%, but usually much better)
CONFIGS = [ int( round( max_cost/float(i) ) ) for i in COSTS ]

for i,params in enumerate(PARAMS):
  print  params, COSTS[i], CONFIGS[i], COSTS[i]*CONFIGS[i]
  time = params[0]*params[1]

  dir="/home/endres/LowEnergyFermions/nrCPS_v2_3_0/production"
  subdir=str(params[0]) + "_" + str(params[1]) + "_" + str(params[2])

  # copy directory
  os.system( "cp -r " +dir + "/two_body-branched " +subdir )

  # replace do.arg time
  os.system( "sed -i \"1c\\"+str(time)+"\" "+subdir+"/args/do.arg" )

  # replace evo.arg configs
  os.system( "sed -i \"2c\\"+str( CONFIGS[i] )+"\" "+subdir+"/args/evo.arg" )

  # replace levels, steps and branches (this is dangerous!!!)
  os.system( "sed -i \"221c\\      int levels = "+str( params[0] )+";\" "+subdir+"/main.C" )
  os.system( "sed -i \"222c\\      int steps = "+str( params[1] )+";\" "+subdir+"/main.C" )
  os.system( "sed -i \"223c\\      int branches = "+str( params[2] )+";\" "+subdir+"/main.C" )

# copy create_binary
os.system( "cp " +dir + "/create_binary.sh ." )

