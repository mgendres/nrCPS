#!/bin/bash

rm -rf callgrind.*
../create_binary.sh
valgrind --tool=callgrind ./a.out 
callgrind_annotate --inclusive=yes callgrind.out.* > callgrind.readable
 ../../scripts/performance/profile.py callgrind.readable
