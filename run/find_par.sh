#!/bin/bash

# submit from project directory !!!


#$ -S /bin/bash
#$ -N ca
#$ -q freya.q
#$ -l h_vmem=6G # job is killed if exceeding this
#$ -cwd
#$ -o ./log/$JOB_NAME.$TASK_ID.$JOB_ID
#$ -j y

export OPENBLAS_NUM_THREADS=1 # avoid multithreading in numpy
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

date

vargs=$(awk "NR==$(($SGE_TASK_ID + 1))" ./run/parameter_list.tsv)
echo "./exe/gh_test_src ${vargs[$id]}"

./exe/gh_test_src ${vargs[$id]}

date
