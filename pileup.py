#!/usr/local/bin/python3
# Plot frequency of shared IBD segments for a set of samples

from collections import Counter, defaultdict
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import sys

def main() :

  if len(sys.argv) != 2 : sys.exit('usage: pileup.py <hmmIBD output (.hmm.txt) file>')
  pdf_file = 'pileup.pdf'
  hmm_file = sys.argv[1]
  binsize = 1000
  
  pp = PdfPages(pdf_file)
  
  # Inefficiently figure out how long the chromsomes are
  chr_ends = Counter()
  max_chrom = 0
  with open(hmm_file, 'r') as inf :
    head = inf.readline().rstrip().split()
    idx = {col: i for i, col in enumerate(head)}
    for line in inf :
      pieces = line.rstrip().split()
      chrom = int(pieces[idx['chr']])
      if chrom > max_chrom : max_chrom += 1
      end = int(pieces[idx['end']])
      if end > chr_ends[chrom] : chr_ends[chrom] = end
    inf.seek(0)

    depth = [ [0]*(chr_ends[x]//binsize+1) for x in range(max_chrom+1) ]

    inf.readline()      # skp the header this time around
    olds1 = olds2 = ''   # previous samples, to tell whether we've started a new one
    ntot_pair = 0
    for line in inf :
      pieces = line.rstrip().split()
      diff = pieces[idx['different']]
      s1 = pieces[idx['sample1']]
      s2 = pieces[idx['sample2']]
      if s1 != olds1 or s2 != olds2 :
        ntot_pair += 1
      olds1 = s1
      olds2 = s2
      if diff == '1' : continue
      chrom = int(pieces[idx['chr']])
      this_clust = None 
      end = int(pieces[idx['end']])
      start = int(pieces[idx['start']])
      bin_end = end // binsize
      bin_start = start // binsize
      for i in range(bin_start, bin_end+1) :
        depth[chrom][i] += 1

  print('total sample pairs', ntot_pair)

  for chrom in range(1, max_chrom+1) :
    xs = [(x+.5)*binsize/1000000 for x in range(len(depth[chrom])) ]   # bin starts
    xc = [(x+.5)*binsize/1000000 for x in range(len(depth[chrom])) ]   # bin centers
    xe = [(x+1)*binsize/1000000 for x in range(len(depth[chrom])) ]   # bin ends

    fig, ax = plt.subplots()
    init1 = False
    frac_depth = np.array(depth[chrom] , dtype=int) / ntot_pair
    for i in range(len(xc)) :
      ax.plot((xs[i],xe[i]), (frac_depth[i], frac_depth[i]), marker='', lw=1.5)
      if init1 :
        ax.plot((xe[i-1],xs[i]), (frac_depth[i-1], frac_depth[i]), marker='', lw=.2)
      if frac_depth[i] > 0 : init1 = True
    ax.set_xlabel('Position on chromosome ' + str(chrom) + ' (Mb)')
    ax.set_ylabel('Fraction of times shared')
    ax.set_title('Sharing pileup, chromosome ' + str(chrom))
#    ax.set_ylim(top=0.8)
    fig.savefig(pp, format='pdf')
    
  pp.close()

main()
