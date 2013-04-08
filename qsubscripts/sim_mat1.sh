#!/bin/bash
#$ -S /bin/bash
#$ -q ahost.q
source /home/rxzhu/.bashrc
time /software/python.2.7.3/bin/python $VECODE/sim_mat_gen.py single 15bits.fp
