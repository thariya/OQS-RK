#!/bin/bash
#PBS -P ty6
#PBS -q expressbw
#PBS -l walltime=01:00:00
#PBS -l mem=256GB
#PBS -l ncpus=28

## The job will be executed from current working directory instead of home.
#PBS -l wd

module load gcc/6.2.0
g++ -fopenmp -I ../Boost/boost_1_61_0 *.cpp -o code.test 
./code.test 5 3 50 1.0 > output.test
