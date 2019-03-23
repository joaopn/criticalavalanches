#!/bin/bash

# submit from project directory !!!


#$ -S /bin/bash
#$ -N ca_256000
#$ -q rostam.q
#$ -l h_vmem=6G # job is killed if exceeding this
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
num_neur=256000
outname=$(printf "$(pwd)/dat/N%s/%05d.hdf5" $num_neur $SGE_TASK_ID)
vargs=(-N\ $num_neur\ -r\ {0..100}\ -s\ $id\ -o\ $outname\ -h\ 1.1e-4\ -m\ 0.9\ -T\ 5e5)
echo "./exe/gh_test_src ${vargs[$id]}"

./exe/gh_test_src ${vargs[$id]}

date
