f = open("data/energy.values", "r")
tpls = map(lambda l: l.split(), f.readlines())
tpls = map(lambda tpl: (tpl[0], float(tpl[1])), tpls)

d = dict(tpls)

cids  = map(lambda s: s[0], tpls)

from lmatrix import lmatrix

m = lmatrix(cids)
for c1 in cids:
    for c2 in cids:
        m[c1,c2] = abs(d[c1] - d[c2])

print m.to_csv_str()        
