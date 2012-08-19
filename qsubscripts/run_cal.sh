#!/bin/bash
#$ -S /bin/bash
#$ -q bhost.q
export SCHRODINGER=/home/rxzhu/software/schrodinger_academic
export SCHRODINGER_PATH=$SCHRODINGER/mmshare-v21025/lib/Linux-x86_64/lib/python2.7/site-packages/
export MMSHARE_EXEC=$SCHRODINGER/mmshare-v21025/bin/Linux-x86_64
export PYTHONPATH=$SCHRODINGER_PATH:~/.local/lib/python2.7/site-packages/
export LD_LIBRARY_PATH=$SCHRODINGER/mmshare-v21025/lib/Linux-x86_64:$LD_LIBRARY_PATH

time /software/python.2.7.3/bin/python /home/rxzhu/code/Virus-Experiment/cal_diff_hydro_vars.py
