# Energy test with optimised permutations

[![DOI](https://zenodo.org/badge/113197119.svg)](https://zenodo.org/badge/latestdoi/113197119)

The [energy test](https://arxiv.org/abs/math/0309164) is a method for determining if two samples orginate from the same underlying distribution. This implimention can efficiently make use of multiple CPUs, with support for using scaled permutations as described in [JINST 13 (2018) no.04, P04011](https://inspirehep.net/record/1648453).

For an implimentation that can be used with NVIDIA GPUs, which may be faster for unscaled calculations with large samples (> 10⁷ points), see [Manet](https://manet.hepforge.org/).

## Compiling

### macOS with gcc

This assumes gcc has been installed using homebrew.

```bash
g++-7 -o energy_test -std=c++1z -O3 -march=native -Wall -I. energy_test.cpp -fopenmp
```

### Manchester machines

```bash
source /cvmfs/sft.cern.ch/lcg/external/gcc/6.2.0/x86_64-slc6-gcc62-opt/setup.sh
g++ -o energy_test -std=c++1z -O3 -march=native -Wall -I. energy_test.cpp -fopenmp
```

### lxplus

```bash
source /cvmfs/sft.cern.ch/lcg/external/gcc/6.2.0/x86_64-slc6-gcc62-opt/setup.sh
g++ -o energy_test -std=c++1z -O3 -Wall -I. energy_test.cpp -fopenmp
```

## Input data

The data is assumed to be in space separated text files and a script is included can create a suitable file (assuming `root_pandas` in installed):

```bash
./prepare_data.py  --input-fn=data/raw/sample0.root --output-fn-D0=data/csv/sample0-d0.csv --output-fn-D0bar=data/csv/sample0-d0bar.csv
```

## Running

Run the energy test for the input data only:

```bash
time ./energy_test sample1.csv sample2.csv
```

Run the energy test for the input data only, limiting the input files to 10000 events:

```bash
time ./energy_test --max-events=10000 sample1.csv sample2.csv
```

Run the energy test with permutations for 10000 points from sample 1 and 11000 points from sample 2:

```bash
time ./energy_test --max-events-1=10000 --max-events-2=11000 --n-permutations=100 sample1.csv sample2.csv
```


# Further Analysis 

To understand and improve upon the energy test, further tests are done. Fundamentally, the energy test consists of a few steps: finding the distances between all the points in each sample as well as across both samples, and then collapsing th sum of these distances into a test statistic T. It is unknown where the performance comes from specifically. 

The goal is to test if we can achieve a comparable level of sensitivity at each step of the process. If so, we have a good idea of where the power of the energy test originates. 

## Runs A

This is fundamentally concerned with the simple distance distributions of each sample; it introduces a few ways of comparing the distance distributions of the two samples (no cross terms - think of the first two terms in the energy test) to see if any statistical significance can be found. The distance we are concerned with is the euclidean distance in the dalitz space - we do not yet use a metric function. 

## Runs B 

This is more gearing towards the energy test - it compares the distance distributions of the first two terms against the distribution of the third term.

## Runs C 

This is essentially applying the energy test, but instead of calculating a single T value for each permutation, we find the T value for each bin, of which each permutation has a certain number of. 
