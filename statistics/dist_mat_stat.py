from nltk import FreqDist
from sim_mat import load_sim_mat
import math
import os
import pylab

from util.load_hydro_vars import load_hydro_var

class CustomFreqDist(FreqDist):
    def __init__(self,samples = None):
        FreqDist.__init__(self,samples)

        self.original_keys = sorted(set(samples))

    def plot_squentially(self, *args, **kwargs):
        if len(args) == 0:
            args = [len(self)]
        samples = self.original_keys
        freqs = [self[sample] for sample in samples]
        ylabel = "Counts"
        # percents = [f * 100 for f in freqs]  only in ProbDist?
        
        pylab.grid(True, color="silver")
        if not "linewidth" in kwargs:
            kwargs["linewidth"] = 2

        if "title" in kwargs: pylab.title(kwargs["title"])
        del kwargs["title"]

        pylab.plot(freqs, **kwargs)
        pylab.xticks(range(len(samples)), [str(s) for s in samples], rotation=90)
        pylab.xlabel("Samples")
        pylab.ylabel(ylabel)
        pylab.show()
    def save_fig(self, *args, **kwargs):
        if len(args) == 0:
            args = [len(self)]
        samples = self.original_keys
        freqs = [self[sample] for sample in samples]
        ylabel = "Counts"
        # percents = [f * 100 for f in freqs]  only in ProbDist?
        
        pylab.grid(True, color="silver")
        if not "linewidth" in kwargs:
            kwargs["linewidth"] = 2

        title = kwargs["title"]
        pylab.title(title)
        del kwargs["title"]

        pylab.plot(freqs, **kwargs)
        pylab.xticks(range(len(samples)), [str(s) for s in samples], rotation=90)
        pylab.xlabel("Samples")
        pylab.ylabel(ylabel)
        
        #F = pylab.gcf()
        #DefaultSize = F.get_size_inches()
        #F.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]*2) )

        pylab.savefig(os.path.join("/home/xiaohan/Desktop/distribution_step_0.02" , title+".png"))

        pylab.hold(False)

def get_freq_dist(matrix , step=0.05):
    return CustomFreqDist([math.floor(c/step) * step for r in matrix for c in r])
            


if __name__ == "__main__":
    mat_ids = load_hydro_var().keys()
    for mat_id in mat_ids:
        mat_id = "%s_dist_mat" %mat_id
        print mat_id
        #mat_id = "WILM950103_dist_mat" 
        sim_mat = load_sim_mat(mat_id)
        dist_mat = sim_mat / sim_mat.diagonal()
        fd = get_freq_dist(dist_mat,step = 0.02)
        fd.save_fig(title = mat_id)
    #fd.tabulate([0.0,0.25,0.5,0.75,1.0])
    #pylab.savefig(os.path.join("/home/xiaohan/Desktop/distribution" , "together.png"))
