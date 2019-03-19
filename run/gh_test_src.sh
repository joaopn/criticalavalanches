#!/bin/bash

#$ -S /bin/bash
#$ -N gh_test_src
#$ -q rostam.q
#$ -cwd
#$ -o /scratch/nst/projects/paul_criticalavalanches/gh_test_src/log/$JOB_NAME.$TASK_ID.$JOB_ID
#$ -j y


export OPENBLAS_NUM_THREADS=1 # avoid multithreading in numpy
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1


# -N for particle numbers
# -r is a dummy arguments for repetitions
# -s seed set to taskid
vargs=(-N\ {144000,72000,36000,18000}\ -r\ {0..50}\ -s\ $SGE_TASK_ID\ -o\ /scratch/nst/projects/paul_criticalavalanches/dat/)
echo ${vargs[$SGE_TASK_ID]}

/scratch/nst/projects/paul_criticalavalanches/run/gh_test_src.py $SGE_TASK_ID
