from time import time,strftime
import os
class Logger(object):
    def __init__(self,name):
        self.name = name

    def log(self,msg):
        raise NotImplementedError("please do something")

class TaskFileLogger(Logger):
    def __init__(self,name):
        Logger.__init__(self , name)
        if not os.path.exists("/home/rxzhu/code/Virus-Experiment/log"):
            os.mkdir("log")
        self.f = open(os.path.join("/home/rxzhu/code/Virus-Experiment/log" , "%s_%s.log" %(self.name,self.get_cur_date_string())),'a')

    def log(self,msg = ""):
        self.f.write("%s : %s at %s\n" %(self.name , msg,self.get_cur_time_string()))
        self.f.flush()

    def get_cur_date_string(self):
        return strftime("%Y-%m-%d")

    def get_cur_time_string(self):
        return strftime("%Y-%m-%d %H:%M:%S")

    def close(self):
        slef.f.close()


