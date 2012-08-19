"""
using the lifehpc.tongji.edu.cn computing nodes to perform calculation tasks
module including task assigment manager
"""

import os
import time
import random

from cal_diff_hydro_vars import do_gen_mat_task,load_hydro_var
from prog_run_states import HydroVarProgStates
from config import *

def fetch_and_do(worker = 'dummy'):
    while True:
        task = fetch_unstarted_task()
        if task is None:
            break
        register_doing_task(task,worker)
        do_gen_mat_task(task)
        #dummy_fun = lambda :time.sleep(1.5)
        #dummy_fun()

        notify_task_done(task,worker)

def init_task_states(task_names = []):
    global db
    col = db["prog_states"]

    col.ensure_index([("key",pymongo.ASCENDING)],unique = True)
    col.save({"key":"unstarted_tasks","value":task_names})

def fetch_unstarted_task():
    global db
    col = db["prog_states"]
    unstarted = col.find_one({"key":"unstarted_tasks"})
            
    col.update({"key":"unstarted_tasks"},{"$pop":{"value":1}})
    try:
        return unstarted["value"][-1]
    except IndexError:#no task left
        print "no task left"
        return None
    
def register_doing_task(task_name,node):
    global db
    col = db["prog_states"]
    if not col.find_one({"key":"doing_tasks"}):
        col.save({"key":"doing_tasks","value":[]})

    col.update({"key":"doing_tasks"},{"$push":{"value":task_name}})
    print "%s will do %s"  %(node,task_name)
    print "unstarted tasks",col.find_one({"key":"unstarted_tasks"})["value"]
    print "doing tasks",col.find_one({"key":"doing_tasks"})["value"]

    db["log"].save({"node":node,"type":"started","time":time.time(),"task_name":task_name})
    
def notify_task_done(task_name,node):
    global db
    col = db["prog_states"]
    col.update({"key":"doing_tasks"},{"$pull":{"value":task_name}})
    if not col.find_one({"key":"done_tasks"}):
        col.save({"key":"done_tasks","value":[]})
        
    col.update({"key":"done_tasks"},{"$push":{"value":task_name}})
    print "done tasks",col.find_one({"key":"done_tasks"})["value"]
    print

    db["log"].save({"node":node,"type":"finished","task_name":task_name,"time":time.time()})

def gen_shell_scripts(script_name  , nodes = []):
    script_tmpl = """
#!/bin/bash
#$ -S /bin/bash
#$ -q %s.q
export SCHRODINGER=/home/rxzhu/software/schrodinger_academic
export SCHRODINGER_PATH=$SCHRODINGER/mmshare-v21025/lib/Linux-x86_64/lib/python2.7/site-packages/
export MMSHARE_EXEC=$SCHRODINGER/mmshare-v21025/bin/Linux-x86_64
export PYTHONPATH=$SCHRODINGER_PATH:~/.local/lib/python2.7/site-packages/
export LD_LIBRARY_PATH=$SCHRODINGER/mmshare-v21025/lib/Linux-x86_64:$LD_LIBRARY_PATH

time /software/python.2.7.3/bin/python /home/rxzhu/code/Virus-Experiment/%s.py -n %s
    """ 

    script_fps = []
    for node in nodes:
        fp = os.path.join(qsub_scripts_path ,"%s_at_%s.sh" %(script_name , node))
        with open(fp,'w') as f:
            f.write(script_tmpl %(node , script_name , node))
        script_fps.append(fp)            
    return script_fps 

def submit_task(script_path):
    os.system("qsub %s" %script_path)

def batch_submit(paths):
    for p in paths:
        time.sleep(1)
        print "dispatching"
        submit_task(p)

if __name__ == "__main__":
    hs = HydroVarProgStates()
    tasks = hs.get_unfinished()

    init_task_states(tasks)

    nodes = ['bhost','ahost','192G']

    fps = gen_shell_scripts(nodes)

    batch_submit(fps)
