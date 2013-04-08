#!/bin/bash
#$ -S /bin/bash
#$ -q ahost.q
source /home/rxzhu/.bashrc
/software/python.2.7.3/bin/python $VE_CODE/sim_mat_gen.py single atg.fp
