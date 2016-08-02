#!/home/endres/Python-2.5.2/python
import os, sys 
sys.path.append("/home/endres/PythonUtils")
from ToFloat import ToFloat 
from Mean import Mean
from Variance import Variance
from Export import Export
from ToString import ToString
from Import import Import
from Power import Power
from Divide import Divide
from Multiply import Multiply
from math import log

# Check some things before processing data
if len(sys.argv)!=3:
  print "Usage: split_data.py <source_stem> <target_stem>"
  sys.exit(-1)
SOURCE_STEM = sys.argv[1]
TARGET_STEM = sys.argv[2]

data = Import(SOURCE_STEM)
data = ToFloat(data)

tol = 0.1
crop = 20

count = 0
split = [data[0]]
for i in range(1,len(data)):
  if (data[i][1]/data[i-1][1] < 0.0) and not (data[i][1]<tol):
    Export(ToString(split), TARGET_STEM+"."+str(count))
    count += 1
    split = []
  else:
    if (abs(data[i][1])<crop):
      split.append(data[i])
Export(ToString(split), TARGET_STEM+"."+str(count))
