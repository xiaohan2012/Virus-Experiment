#!/usr/bin/env python

import os,sys,optparse

def prepend_zeros(tab, aa_diff):
  return [0 for x in range(aa_diff*9)] + tab
  


def gen_avg_sift(inname,oname):
    cutoff = 0.5
    avg_fp = []
    min_aa = 0
    fp_count = 0
    
    sift_fh = open(inname)
    for line in sift_fh:
        rec, lig, aa, fp = line.strip().split(':')
        if min_aa == 0: #reading first fp
          min_aa = int(aa)
          avg_fp = list(fp)
          fp_count += 1
          continue
        if int(aa) < min_aa:
          avg_fp = prepend_zeros(avg_fp, min_aa-int(aa))
          min_aa = int(aa)
        if int(aa) > min_aa:
          fp = prepend_zeros(list(fp), int(aa)-min_aa)
        avg_fp = [(int(x) + int(y)) for x, y in zip(avg_fp, fp)]
        fp_count += 1
    sift_fh.close()

    avg_fp = map(lambda x: round(float(x)/fp_count,2) if float(x)/fp_count >= cutoff else 0, avg_fp)


    pattern_fh = open(oname,'w')
    pattern_fh.write(str(min_aa)+ ':'+' '.join(str(x) for x in avg_fp) + "\n\n")
    for chunk_no in range(len(avg_fp[::9])):
      pattern_fh.write("%i\t\t\t%s\n" % (int(chunk_no)+int(min_aa), ' '.join(str(x) for x in avg_fp[chunk_no*9:(chunk_no+1)*9])))

    print 'write output pattern to %s' %oname
    pattern_fh.close()
