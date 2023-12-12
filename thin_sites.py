#!/usr/bin/env python3

from collections import Counter
import sys

def main() :
  min_space = 50
  min_freq = 0.005
  binsize = 2000
  max_per_bin = 12
  if len(sys.argv) != 3 : sys.exit('usage: thin_sites.py <input freq file> <output file of good sites>')
  outf = open(sys.argv[2], 'w')
  
  max_freq_lists = {x : [] for x in range(15)}  # the largest minor allele frequency for this snp (=minor allele frequency for biallelic)
  pos_lists = {x : [] for x in range(15)}
  chr_ends = Counter()
  with open(sys.argv[1], 'r') as inf :
    for line in inf :
      pieces = line.rstrip().split()
      chrom = int(pieces[0])
      pos = int(pieces[1])
      if pos > chr_ends[chrom] : chr_ends[chrom] = pos
      pos_lists[chrom].append(pos)
      maxf = max([float(x) for x in pieces[2:]])
      if maxf > 0.5 : maxf = 1 - maxf
      max_freq_lists[chrom].append( maxf )

  count_sum = 0
  for chrom in range(1,15) :
    snp_count = [0]*(chr_ends[chrom]//binsize+1)
    chr_count = 0     # kept on this chrom
    max_freqs = max_freq_lists[chrom]
    poss = pos_lists[chrom]
    nsnp = len(max_freqs)
    to_use = [False] * nsnp
    sorted_index = sorted(range(nsnp), key=lambda k: max_freqs[k], reverse=True)
    for idx in sorted_index :
      if max_freqs[idx] < min_freq : break
      posi = poss[idx]
      ibin = posi // binsize
      if snp_count[ibin] >= max_per_bin :
        continue
      # First look on left side
      skip_it = False
      for jdx in range(idx-1, -1, -1) :
        posj = poss[jdx]
        if posi - posj > min_space : break
        if to_use[jdx] :    # this neighboring SNP is already in use, so I'll have to drop idx
          skip_it = True
          break
      if skip_it : continue
      # Then right side
      for jdx in range(idx+1, len(poss)) :
        posj = poss[jdx]
        if posj - posi > min_space : break
        if to_use[jdx] :    # this neighboring SNP is already in use, so I'll have to drop idx
          skip_it = True
          break
      if skip_it : continue
      snp_count[ibin] += 1
      to_use[idx] = True
      chr_count += 1
    print(chrom, 'N snp:', nsnp, 'kept:', chr_count)
    count_sum += chr_count
    for kdx in range(nsnp) :
      if to_use[kdx] : print('{:d}\t{:d}'.format(chrom, poss[kdx]), file=outf)
  print(count_sum)
    
main()
