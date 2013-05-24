#!/bin/bash
#$ -S /bin/bash
#$ -q bhost.q
source /home/rxzhu/.bashrc
source /home/rxzhu/.bashrc_xiaohan
time /software/python.2.7.3/bin/python $VE_CODE/fp/fp_16.py comp_plane
