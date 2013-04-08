#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N sim
#$ -pe mpi 1
#$ -j y
#$ -q ahost.q


source /home/rxzhu/.bashrc

/software/python.2.7.3/bin/python $VE_CODE/sim_mat_gen.py single atb.fp 
