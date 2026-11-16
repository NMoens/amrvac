#!/bin/bash
cd /home/nicolasm/Codes/amrvac3/amrvac/tests/hd/m1_Ostar_1D
for n in Euler_tvdlf_minmod Midpoint_tvdlf_minmod Euler_hll_minmod Euler_tvdlf_koren Euler_tvdlf_weno3; do
  (mpirun -n 3 ./amrvac -i solv/$n.par > solv/$n.out 2>&1; echo "$n exit $? $(date +%T)" >> solv/status.txt) &
done
wait
echo ALLDONE >> solv/status.txt
