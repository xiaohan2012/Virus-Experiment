#!/bin/bash
#$ -S /bin/bash
#$ -q bhost.q
source /home/rxzhu/.bashrc
source /home/rxzhu/.bashrc_xiaohan

time /software/python.2.7.3/bin/python $VE_CODE/fp/fp_370.py atb
/software/python.2.7.3/bin/python $VE_CODE/qsubscripts/notifier.py "fp_370 atb as rec"

