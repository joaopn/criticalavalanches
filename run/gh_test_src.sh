#!/bin/bash

# submit from project directory !!!

#$ -S /bin/bash
#$ -N gh_test_src
#$ -q freya.q
#$ -cwd
#$ -o ./log/$JOB_NAME.$TASK_ID.$JOB_ID
#$ -j y


export OPENBLAS_NUM_THREADS=1 # avoid multithreading in numpy
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

date

# -N for particle numbers
# -r is a dummy arguments for repetitions
# -s seed set to taskid
id=$(($SGE_TASK_ID - 1))
vargs=(-N\ {144000,72000,36000,18000}\ -r\ {0..50}\ -s\ $id\ -o\ ./dat/\ -h\ 1.1e-4\ -m\ 0.9\ -T\ 1e5)
echo $(./exe/gh_test_src $(echo ${vargs[$id]}))

./exe/gh_test_src $(echo ${vargs[$id]})

date
