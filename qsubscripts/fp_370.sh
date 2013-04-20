#!/bin/bash
#$ -S /bin/bash
#$ -q bhost.q
source /home/rxzhu/.bashrc
time /software/python.2.7.3/bin/python $VE_CODE/fp/fp_370.py
