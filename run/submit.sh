#!/bin/bash

# submit from project directory !!!

#$ -S /bin/bash
#$ -N tpl
#$ -q rostam.q
#$ -l h_vmem=6G # job is killed if exceeding this
#$ -cwd
#$ -o ./log/
#$ -j y

export OPENBLAS_NUM_THREADS=1 # avoid multithreading in numpy
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

# notify "started $JOB_NAME.$TASK_ID.$JOB_ID"
date

vargs=$(awk "NR==$(($SGE_TASK_ID + 1))" ./run/parameters.tsv)
echo "${vargs[$id]}"

${vargs[$id]}

date
# notify "finished $JOB_NAME.$TASK_ID.$JOB_ID"
