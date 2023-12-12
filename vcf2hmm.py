#!/usr/bin/env python3

# Script to output genotype data from VCF file formatted for hmmIBD

from collections import Counter
import gzip
import argparse
import re
import sys

def main() :
  # Accept only sites with filter state = 'PASS'
  filter_sites = False
  kill_indel = False
  # Mininum genotyping call rate to accept variant
  min_call = 0.80
  min_depth = 0      # minimum read depth to accept a genotype
  # Option: minimum number of minor allele copies to keep variant
  #  (1 = kill monomorphic, 2 = kill singletons)
  min_copy = 0
  parser = argparse.ArgumentParser(description='Extract genotypes from VCF file for input to hmmIBD')
  parser.add_argument('vcf_file', help='input file name (.vcf or .vcf.gz)')
  parser.add_argument('out_file', help='text file name base (omit extensions)')
  parser.add_argument('-s', '--samp_file', help='Name of file containing names of samples to keep')
  parser.add_argument('-l', '--loci_file', help='Name of file containing variant loci (chrom pos) to keep')
  args = parser.parse_args()
  infile = args.vcf_file
  filebase = args.out_file
  seqfile = f'{filebase}_seq.txt'
  freqfile = f'{filebase}_freq.txt'
  allfile = f'{filebase}_allele.txt'
  sampfile = args.samp_file
  snpfile = args.loci_file
  
  all_samps = True
  sstr = 'all'
  good_samps = set()
  if sampfile != None :
    all_samps = False
    sstr = sampfile
    with open(sampfile, 'r') as sampf :
      for line in sampf :
        samp = line.rstrip()
        good_samps.add(samp)
  print('Sample source:', sstr)
  
  all_snps = True
  sstr = 'all'
  good_snps = set()
  if snpfile != None :
    all_snps = False
    sstr = snpfile
    with open(snpfile, 'r') as snpf :
      for line in snpf :
        pieces = line.rstrip().split()
        chrom = int(pieces[0])
        pos = int(pieces[1])
        good_snps.add( (chrom, pos) )
  print('SNP source:', sstr)

  # Decide whether the file is compressed or not, and open it accordingly
  if re.search(r'\.gz$', infile) :
    inf = gzip.open(infile, 'rt')
  else :
    inf = open(infile, 'r')

  samples = []
  nindel = n_lowcall = n_lowfreq = 0
  n_fail_filter = n_pass_filter = 0
  chrom_name_map = {}
  n_chrom_names = 0
  nline = 0
  with open(seqfile, 'w') as seqf, open(allfile, 'w') as allf, open(freqfile, 'w') as freqf :
    for line in inf :
      if re.match(r'\#CHROM', line) :
        samples = line.rstrip().split('\t')[9:]
        print(f'samples: {len(samples)}, of which ', end='')
        for samp in samples :
          if all_samps :
            subp = samp.split('_')
            if len(subp) < 3 or subp[0] != 'SEN' :
              pass
              # continue      # non-standard Senegal name
            if samp in good_samps : print('Duplicate sample', samp)
          elif samp not in good_samps :
            continue
          good_samps.add(samp)
        print(f'{len(good_samps)} were kept')
        print('chrom\tpos', sep='', end='', file=seqf)
        for samp in samples :
          if samp not in good_samps : continue
          print('\t', samp, sep='', end='', file=seqf)
        print('', file=seqf)
        continue
      elif re.match(r'\#', line) : continue
      
      # data line
      nline += 1
      pieces = line.rstrip().split('\t')
      chrom_str = pieces[0]
      filter_state = pieces[6]
      # Convert chrom names into integers
      if chrom_str in chrom_name_map :
        chrom = chrom_name_map[chrom_str]
      else :
        match = re.search(r'^Pf3D7_(\d+)_v3$', chrom_str)
        if match :
          chrom = int(match.group(1))
          chrom_name_map[chrom_str] = chrom
          print(f'chromosome {chrom_str} maps to {chrom}')
          n_chrom_names += 1
        else :
          match2 = re.search(r'^\d{1,2}$', chrom_str, re.A)
          if match2 :
            chrom = int(chrom_str)
            chrom_name_map[chrom_str] = chrom
            print(f'chromosome {chrom_str} maps to {chrom}')
            n_chrom_names += 1
          else :
            n_chrom_names += 1
            chrom = n_chrom_names
            chrom_name_map[chrom_str] = chrom
            print(f'chromosome {chrom_str} maps to {chrom}')

      pos = int(pieces[1])
      if not all_snps and (chrom, pos) not in good_snps : continue
      ref_all = pieces[3]
      alt_alls = pieces[4].split(',')
      kill = False
      if not re.match(r'^[ACGT\.]$', ref_all) :
        if kill_indel :
          nindel += 1
          continue
      for alt_all in alt_alls :
        if not re.match(r'[ACGT\.]$', alt_all) :
          kill = True
      if kill :
        if kill_indel :
          nindel += 1
          continue
      if filter_sites and filter_state != 'PASS' :
        n_fail_filter += 1
        continue
      n_pass_filter += 1

      formats = pieces[8].split(':')
      genotypes = pieces[9:]
      outline = '{:d}\t{:d}'.format(chrom, pos)
      nassay = ncall = 0

      all_calls = Counter()
      for isamp, samp in enumerate(samples) :
        if samp not in good_samps : continue
        genotype = genotypes[isamp].split(':')
        call = ''
        depth = 0
        for iform, form in enumerate(formats) :
          if form == 'GT' :
            call = genotype[iform]
          elif form == 'DP' :
            if genotype[iform] != '.' : 
              depth = int(genotype[iform])
        allele = '-2'
        if re.match(r'[0-9]+', call) and depth >= min_depth :
          match = re.match(r'([0-9]+)[\/|]([0-9]+)', call)
          g1 = match.group(1)
          g2 = match.group(2)
          if g1 != g2 :
            # het
            allele = '-1'
          else :
            allele = g1
            all_calls[int(allele)] += 1
            ncall += 1
        else :
          #nocall
          allele = '-1'
        nassay += 1
        outline += ('\t' + allele)
      if nassay == 0 or ncall / nassay < min_call :
        n_lowcall += 1
        continue
      max_count = max(list(all_calls.values()))
      max_allele = max(list(all_calls.keys()))
      if ncall - max_count < min_copy :
        n_lowfreq += 1
        continue
      print(outline, file=seqf)
      print('{:d}\t{:d}'.format(chrom, pos), ref_all, '\t'.join(alt_alls), sep='\t', file=allf)
      print('{:d}\t{:d}'.format(chrom, pos), end='', file=freqf)
      for iall in range(max_allele+1) :
        print(f'\t{all_calls[iall] / ncall : .4f}', end='', file=freqf)
      print('', file=freqf)
  print('N indels killed', nindel)
  print('N failed/passed filter', n_fail_filter, n_pass_filter)
  print('N dropped for low call rate', n_lowcall)
  print('N dropped for low minor allele freq', n_lowfreq)
main()
