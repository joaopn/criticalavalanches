#!/bin/bash

# submit from project directory !!!

#$ -S /bin/bash
#$ -N gh_test_src
#$ -q freya.q
#$ -l h_vmem=5G
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
# -s seed set to taskid - 1, for some reason qsub doesnt want array_id 0
id=$(($SGE_TASK_ID - 1))
outname=$(printf "./dat/%05d.hdf5" $SGE_TASK_ID)
vargs=(-N\ 256000\ -r\ {0..50}\ -s\ $id\ -o\ $outname\ -h\ 1.1e-4\ -m\ 0.9\ -T\ 5e5)
echo "./exe/gh_test_src ${vargs[$id]}"

./exe/gh_test_src ${vargs[$id]}

date
