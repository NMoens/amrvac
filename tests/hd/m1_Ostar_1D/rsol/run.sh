#!/bin/bash
cd /home/nicolasm/Codes/amrvac3/amrvac/tests/hd/m1_Ostar_1D
echo "start $(date +%T)" > rsol/status.txt
run(){ mpirun -n $2 ./amrvac -i rsol/$1.par > rsol/$1.out 2>&1; echo "$1 exit $? $(date +%T)" >> rsol/status.txt; }
run R0.1 3 & run R0.3 5 & run R1.0 8 &
wait; echo ALLDONE >> rsol/status.txt
