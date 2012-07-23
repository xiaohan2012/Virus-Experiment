h_fp = "hydro_variations.dat"
def load_hydro_var(h_fp):
    with open(h_fp,'r') as f:
        reses = f.readline().split()
        d_ = {}
        for l in f.readlines():
            linkdb = l.split()[0]
            d_[linkdb] = {}
            for res,num in zip(reses,l.split()[1:]):
                d_[linkdb][res] = float(num)
    return d_

if __name__ == "__main__":
    hydro_vars = load_hydro_var(h_fp)
    print hydro_vars 

