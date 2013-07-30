#!/bin/bash
#$ -S /bin/bash
#$ -q bhost.q
source /home/rxzhu/.bashrc
source /home/rxzhu/.bashrc_xiaohan
/software/python.2.7.3/bin/python $VE_CODE/fp/fp_175_padded.py complex atg res 
/software/python.2.7.3/bin/python $VE_CODE/fp/fp_175_padded.py complex atg tri

/software/python.2.7.3/bin/python $VE_CODE/fp/fp_175_padded.py complex atb res 
/software/python.2.7.3/bin/python $VE_CODE/fp/fp_175_padded.py complex atb tri
 
/software/python.2.7.3/bin/python $VE_CODE/fp/fp_175_padded.py residue atg res 
/software/python.2.7.3/bin/python $VE_CODE/fp/fp_175_padded.py residue atg tri 

/software/python.2.7.3/bin/python $VE_CODE/fp/fp_175_padded.py residue atb res 
/software/python.2.7.3/bin/python $VE_CODE/fp/fp_175_padded.py residue atb tri 
