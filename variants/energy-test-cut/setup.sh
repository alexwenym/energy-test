source /cvmfs/sft.cern.ch/lcg/external/gcc/6.2.0/x86_64-slc6-gcc62-opt/setup.sh
g++ -o energy_test -std=c++1z -O3 -Wall -I. energy_test.cpp -fopenmp
