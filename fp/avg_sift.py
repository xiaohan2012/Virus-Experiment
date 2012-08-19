#!/usr/bin/env python
"""
generate the average finger print file from one or more finger-print files
"""
import os,sys,optparse

def prepend_zeros(tab, aa_diff):
  return ["00" for x in range(aa_diff*13)] + tab
  

def get_avg_fp(fp):
    for i in xrange(len(fp)/2):
        i1 = 2*i
        i2 = i1 + 1
        yield fp[i1] + fp[i2]
        
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
          #avg_fp = list(fp)
          avg_fp = list(get_avg_fp(fp))
          #print avg_fp
          fp_count += 1
          continue
        #the two "if"s are redundant
        if int(aa) < min_aa:
          avg_fp = prepend_zeros(avg_fp, min_aa-int(aa))
          min_aa = int(aa)
        if int(aa) > min_aa:
          fp = prepend_zeros(list(fp), int(aa)-min_aa)
        avg_fp = [(int(x) + int(y)) for x, y in zip(avg_fp, fp)]
        fp_count += 1
    sift_fh.close()

    #avg_fp = map(lambda x: round(float(x)/fp_count,2) if float(x)/fp_count >= cutoff else 0, avg_fp)


    pattern_fh = open(oname,'w')
    pattern_fh.write(str(min_aa)+ ':'+' '.join(str(x) for x in avg_fp) + "\n\n")
    for chunk_no in range(len(avg_fp[::13])):
      pattern_fh.write("%i\t\t\t%s\n" % (int(chunk_no)+int(min_aa), ' '.join(str(x) for x in avg_fp[chunk_no*13:(chunk_no+1)*13])))

    #print 'write output pattern to %s' %oname
    pattern_fh.close()
