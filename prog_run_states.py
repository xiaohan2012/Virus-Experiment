import pymongo as pm
from cal_diff_hydro_vars import load_hydro_var

class HydroVarProgStates(object):
    
    def __init__(self):
        self.all_var_names = load_hydro_var().keys()
        self.total = 81003
        self.conn = pm.Connection("10.10.20.200")
        self.db = self.conn["virus_cluster"]
        self.states = {}
        self.all = ["%s_dist_mat" %name for name in self.all_var_names]

        self.refresh()

    def _get_stat(self):
        for c_name , h_name in zip(self.all,self.all_var_names):
            col = self.db[c_name]
            self.states[h_name] = col.find().count()
        print self.states

    def get_all(self):
        return self.all

    def get_unfinished(self):
        self.unfinished = []
        for h_name,f_count in self.states.items():
            if f_count != self.total:
                self.unfinished.append(h_name)
        self.unfinished = sorted(self.unfinished,key = lambda a:self.states[a],reverse = True)
        return self.unfinished

    def get_finished(self):
        self.finished = []
        for h_name,f_count in self.states.items():
            if f_count == self.total:
                self.finished.append(h_name)
        return self.finished                

    def display_progress(self):
        self.refresh()
        
        fin_share = len(self.finished) * self.total + sum([self.states[unfin] for unfin in self.unfinished])
        total_share = ( len(self.finished) + len(self.unfinished) ) * self.total

        print "%d out of %d are finished, finish rate %3.2f%%" %(len(self.finished) , len(self.all) , \
                                                                         fin_share / float(total_share)  * 100)

        for fin in self.finished:
            print "%s progress:%s" %(fin , '#' * 20)

        for unfin in self.unfinished:
            print "%s progress:%s" %(unfin , '#' * int(20. * self.states[unfin] / self.total))
        
        print "\nThanks for your patience, we are doing better and better"

    def refresh(self):
        self._get_stat()
        self.get_finished()
        self.get_unfinished()

if __name__ == "__main__":
    hs = HydroVarProgStates()
    hs.display_progress()
