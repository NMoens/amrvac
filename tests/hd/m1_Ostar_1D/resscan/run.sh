#!/bin/bash
# 1D radial resolution scan for Star8 (see base.par). Runs all resolutions in
# parallel from tests/hd/m1_Ostar_1D; progress in resscan/status.txt,
# stdout per run in resscan/nXXXX.out.
#   bash resscan/run.sh            (all resolutions)
#   bash resscan/run.sh 0640 1280  (a subset)
set -u
cd "$(dirname "$0")/.." || exit 1
make -j8 > resscan/make.log 2>&1 || { echo "build failed, see resscan/make.log"; exit 1; }

# MPI ranks per resolution (blocks of 40 cells: 8, 16, 32, 64 blocks)
declare -A NP=( [0320]=2 [0640]=2 [1280]=4 [2560]=8 )
runs=("${@:-0320 0640 1280 2560}")
runs=(${runs[@]})

echo "start $(date +%F_%T)" >> resscan/status.txt
for n in "${runs[@]}"; do
  (
    mpirun -n "${NP[$n]}" ./amrvac -i resscan/base.par "resscan/n$n.par" > "resscan/n$n.out" 2>&1
    echo "n$n exit $? $(date +%F_%T)" >> resscan/status.txt
  ) &
done
wait
echo "all done $(date +%F_%T)" >> resscan/status.txt
