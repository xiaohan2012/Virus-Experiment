#!/bin/bash
#$ -S /bin/bash
#$ -q bhost.q
source /home/rxzhu/.bashrc

time /software/python.2.7.3/bin/python /home/rxzhu/QiuTianyi/code/ve/fp/epi_gen.py 1> /home/rxzhu/QiuTianyi/log/paraepi_gen.o.log 2> /home/rxzhu/QiuTianyi/log/paraepi_gen.e.log
