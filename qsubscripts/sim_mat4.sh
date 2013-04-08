#!/bin/bash
#$ -S /bin/bash
#$ -q ahost.q
source /home/rxzhu/.bashrc
time /software/python.2.7.3/bin/python $VE_CODE/sim_mat_gen.py double 15bits.fp
