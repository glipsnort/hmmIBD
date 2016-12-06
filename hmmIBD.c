// Calculate segments of haploid genome identity between pairs of samples. 
// Can filter out sample pairs with high discordance (based on sites where at least 
// 1 of the genomes has the minor allele) to skip pairs with little or no sharing
// States: 0=IBD (identical by descent), 1=non-identical
//  HMM notation follows Rabiner (Proc of the IEEE, V77, p257 (1989)).
//  Note: alpha/beta underflow can be handled either by working in log space + logsumexp trick, or 
//   by scaling (as in Rabiner). I use scaling.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <unistd.h>

int main(int argc, char **argv) {
  /* User-settable parameters */
  const double eps = .001;          // error rate in genotype calls
  const int min_inform = 0;       // minimum number of informative sites (those w/ minor allele)
  const double min_discord = 0.0; // minimum discordance; set > 0 to skip identical pairs
  const double max_discord = 1.0;  // set < 1 to skip unrelated pairs
  const int nchrom = 14;            // 14 for falciparum
  const int min_snp_sep = 10;       // skip next snp(s) if too close to last one; in bp
  const double rec_rate = 5.8e-7;   // recombination rate:5.8e-5 cM/bp, or 17kb/cM
  const double fit_thresh_dpi = .001;
  const double fit_thresh_dk = .1;
  const double fit_thresh_drelk = .001;
  const double k_rec_init = 1.0;          // starting value for N generation parameter

  /* end user-settable parameters */
  double k_rec;  // recombination coefficient (= effective number of meioses)
  const int max_all = 8;
  int niter = 5;    // maximum number of iterations of fit; can be overridden with option -n
  int max_snp = 30000;
  char data_file[128];
  char out_filebase[128], freq_file[128], good_file[128], bad_file[128];
  int max_bad = 100;
  int max_good = 200;
  int linesize = 4000;
  char *newLine, *token, *running, **sample, *head;
  char file[64], **bad_samp=NULL, **good_pair[2]={NULL};
  int itoken, nsample=0, isamp, chr, sum, iall, maxall, all, js, snp_ind, **geno;
  double **discord, pright, seq_ibd_fb=0, seq_dbd_fb=0, p_ibd;
  double **freq=NULL, *ffreq=NULL, xisum, xi[2][2], trans_pred, trans_obs;
  double *phi[2], pinit[2], pi[2], *b[2], a[2][2], ptrans, *alpha[2], *beta[2], *scale;
  double maxval, max_phi=0, max_phiL, seq_ibd, seq_dbd, count_ibd_fb, count_dbd_fb;
  double gamma[2], last_pi=0, last_prob=0, last_krec=0, delpi, delprob, delk;
  FILE *inf=NULL, *outf=NULL, *pf=NULL, *ff=NULL;
  int *diff=NULL, *same_min=NULL, jsamp, *allcount, *use_sample, *traj, add_seq;
  int nsample_use, nsnp, ipair, npair, isnp, chrlen, *pos, *psi[2], max, iline;
  int *nmiss_bypair=NULL, totall, *start_chr=NULL, *end_chr=NULL, is, maxlen;
  int **use_pair=NULL, *nall=NULL, killit, nuse_pair=0, gi, gj, delpos;
  int ntri=0, ibad, nbad, start_pos, ex_all=0, majall, last_snp, c, iflag, oflag;
  int freq_flag, fpos=0, fchr=0, iter, ntrans, nflag, finish_fit, bflag, gflag;
  int prev_chrom, ngood;

  pinit[0] = 0.5;  // flat prior
  pinit[1] = 0.5;

  char usage_string[512];
  strcpy(usage_string, "Usage: hmmIBD -i <input filename> -o <output filename>");
  strcat(usage_string, " [-n <max fit iteration>]\n");
  strcat(usage_string, "     [-f <allele frequency file>] [-b <file with samples to skip>]");
  strcat(usage_string, "  [-g <file with sample pairs to use>]\n");

  opterr = 0;
  nflag = iflag = oflag = freq_flag = bflag = gflag = 0;
  while ( (c = getopt(argc, argv, ":f:i:o:n:b:g:")) != -1) {
    switch(c) {
    case 'f':
      freq_flag = 1;
      strcpy(freq_file, optarg);
      break;
    case 'b':
      bflag = 1;
      strcpy(bad_file, optarg);
      break;
    case 'g':
      gflag = 1;
      strcpy(good_file, optarg);
      break;
    case 'n':
      nflag = 1;
      niter = strtol(optarg, NULL, 10);
      break;
    case 'i':
      iflag = 1;
      strcpy(data_file, optarg);
      break;
    case 'o':
      oflag = 1;
      strcpy(out_filebase, optarg);
      break;
    case ':':
      fprintf(stderr, "option %c requires an argument\n", optopt);
      exit(0);
    case '?':
      fprintf(stderr, "%s",usage_string);
      exit(0);
    }
  }

  if (optind != argc || iflag == 0 || oflag == 0) {
    fprintf(stderr, "%s", usage_string);
    exit(0);
  }

  if (freq_flag == 1) {
    ff = fopen(freq_file, "r");
    if (ff == NULL) {fprintf(stderr, "Could not open frequency file %s\n", freq_file); exit(0);}
  }

  allcount = malloc((max_all+1) * sizeof(int));
  bad_samp = malloc(max_bad * sizeof(char*));
  good_pair[0] = malloc(max_good * sizeof(char*));
  good_pair[1] = malloc(max_good * sizeof(char*));
  newLine = malloc((linesize+1) * sizeof(char));
  assert(newLine != NULL);
  nall = malloc(max_snp * sizeof(int));
  freq = malloc(max_snp * sizeof(double*));
  pos = malloc(max_snp * sizeof(int));
  start_chr = malloc((nchrom+1) * sizeof(int));
  end_chr = calloc(nchrom+1, sizeof(int));
  for (isnp = 0; isnp < max_snp; isnp++) {
    freq[isnp] = malloc((max_all+1) * sizeof(double));
  }
  ffreq = malloc((max_all+1) * sizeof(double));
  for (chr = 1; chr <= nchrom; chr++) {start_chr[chr] = 10000000;}
  for (isamp = 0; isamp < max_bad; isamp++) {
    bad_samp[isamp] = malloc(64 * sizeof(char));
    good_pair[0][isamp] = calloc(64, sizeof(char));
    good_pair[1][isamp] = calloc(64, sizeof(char));
  }

  ngood = nbad = 0;
  if (bflag == 1) {
    inf = fopen(bad_file, "r");
    nbad = 0;
    if (inf == NULL) {
      fprintf(stderr, "Could not open file of bad samples: %s\n", bad_file);
      bflag = 0;
    }
    else {
      while (fgets(newLine, linesize, inf) != NULL) {
	newLine[strcspn(newLine, "\r\n")] = 0;  
	strncpy(bad_samp[nbad], newLine, 63);
	nbad++;
	if (nbad == max_bad) {
	  bad_samp = realloc(bad_samp, 2*max_bad*sizeof(char*));
	  for (isamp = max_bad; isamp < 2*max_bad; isamp++) {
	    bad_samp[isamp] = malloc(64 * sizeof(char));
	  }
	  max_bad *= 2;
	}
      }
      fclose(inf);
    }
  }
  if (gflag == 1) {
    inf = fopen(good_file, "r");
    ngood = 0;
    if (inf == NULL) {
      fprintf(stderr, "Could not open file of good pairs: %s\n", good_file);
      gflag = 0;
    }
    else {
      while (fgets(newLine, linesize, inf) != NULL) {
	newLine[strcspn(newLine, "\r\n")] = 0;  
	for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL && itoken < 2; 
	     itoken++) {
	  strncpy(good_pair[itoken][ngood], token, 63);
	}
	ngood++;
	if (ngood == max_good) {
	  good_pair[0] = realloc(good_pair[0], 2*max_good*sizeof(char*));
	  good_pair[1] = realloc(good_pair[1], 2*max_good*sizeof(char*));
	  for (isamp = max_good; isamp < 2*max_good; isamp++) {
	    good_pair[0][isamp] = calloc(64, sizeof(char));
	    good_pair[1][isamp] = calloc(64, sizeof(char));
	  }
	  max_good *= 2;
	}
      }
      fclose(inf);
    }
  }
  
  inf = fopen(data_file, "r");
  if (inf == NULL) {fprintf(stderr, "Could not open input file %s\n", data_file); exit(0);}
  sprintf(file, "%s.hmm.txt", out_filebase);
  outf = fopen(file, "w");
  if (outf == NULL) {fprintf(stderr, "Could not open output file %s\n", file); exit(0);}
  fprintf(outf, "sample1\tsample2\tchr\tstart\tend\tdifferent\n");
  sprintf(file, "%s.hmm_fract.txt", out_filebase);
  pf = fopen(file, "w");
  if (pf == NULL) {fprintf(stderr, "Could not open output file %s\n", file); exit(0);}
  fprintf(pf, "sample1\tsample2\tN_informative_sites\tdiscordance\tlog_p\tN_fit_iteration\tN_generation");
  fprintf(pf, "\tN_state_transition\tseq_shared_best_traj\tfract_sites_IBD\n");
  
  fgets(newLine, linesize, inf); // header
  while (strlen(newLine) > linesize-2) {
    fseek(inf, 0, 0);
    linesize *= 2;
    free(newLine);
    newLine = malloc((linesize+1) * sizeof(char));
    fgets(newLine, linesize, inf); // header
  }
  newLine[strcspn(newLine, "\r\n")] = 0;  
  head = malloc((linesize+1) * sizeof(char));
  assert(head != NULL);
  strcpy(head, newLine);
  for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
    if (itoken > 1) {nsample++;}
  }
  
  sample = malloc(nsample * sizeof(char*));
  assert(sample != NULL);
  geno = malloc(nsample * sizeof(int*));
  assert(geno != NULL);
  use_sample = malloc(nsample * sizeof(int));
  assert(use_sample != NULL);
  discord = malloc(nsample * sizeof(double*));
  assert(discord != NULL);
  use_pair = malloc(nsample * sizeof(int*));
  assert(use_pair != NULL);
  for (isamp = 0; isamp < nsample; isamp++) {
    geno[isamp] = malloc(max_snp * sizeof(int));
    assert(geno[isamp] != NULL);
    sample[isamp] = malloc(64 * sizeof(char));
    assert(sample[isamp] != NULL);
    use_pair[isamp] = calloc(nsample, sizeof(int));
    assert(use_pair[isamp] != NULL);
    discord[isamp] = malloc(nsample * sizeof(double));
    assert(discord[isamp] != NULL);
  }
  
  isamp = nsample_use = 0;
  prev_chrom = -1;
  for (running = head, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
    if (itoken > 1) {
      strncpy(sample[isamp], token, 63);
      use_sample[isamp] = 1;
      for (ibad = 0; ibad < nbad; ibad++) {
	if (strcmp(token, bad_samp[ibad]) == 0) {
	  use_sample[isamp] = 0;
	  fprintf(stdout, "killing sample %s\n", token);
	  break;
	}
      }
      if (use_sample[isamp] == 1) {
	nsample_use++;
      }
      isamp++;
    }
  }

  fprintf(stdout, "Minimum marker spacing (bp): %d\n", min_snp_sep);
  fprintf(stdout, "Minimum informative markers: %d\n", min_inform);
  fprintf(stdout, "Pairs accepted with discordance in range (%.2f%%, %.2f%%)\n", 
	  min_discord*100, max_discord*100);
  fprintf(stdout, "Genotyping error rate: %.2f%%\n", eps*100);
  fprintf(stdout, "Input file: %s\n", data_file);
  fprintf(stdout, "Frequency file: ");
  if (freq_flag == 1) {fprintf(stdout, "%s\n", freq_file);}
  else {fprintf(stdout, "none\n");}
  npair = nsample_use * (nsample_use-1) / 2;
  fprintf(stdout, "nsample: %d\t used: %d expected pairs: %d\n", nsample, nsample_use, 
	  npair);

  nmiss_bypair = calloc(npair, sizeof(int));
  same_min = calloc(npair, sizeof(int));
  diff = calloc(npair, sizeof(int));

  // All pairs are to be used by default, unless we've read a file of good pairs
  for (isamp = 0; isamp < nsample; isamp++) {
    for (jsamp = isamp+1; jsamp < nsample; jsamp++) {
      if (gflag == 1) {
	use_pair[isamp][jsamp] = 0;
	for (ipair = 0; ipair < ngood; ipair++) {
	  if ( (strcmp(good_pair[0][ipair], sample[isamp]) == 0 && 
		strcmp(good_pair[1][ipair], sample[jsamp]) == 0) ||
	       (strcmp(good_pair[0][ipair], sample[jsamp]) == 0 && 
		strcmp(good_pair[1][ipair], sample[isamp]) == 0) ) {
	    use_pair[isamp][jsamp] = 1;
	  }
	}
      }
      else {
	use_pair[isamp][jsamp] = 1;
      }
    }
  }

  nsnp = 0;
  iline = -1;
  while (fgets(newLine, linesize, inf) != NULL) {
    if (nsnp == max_snp) {
      nall = realloc(nall, 2*max_snp*sizeof(int));
      assert(nall != NULL);
      pos = realloc(pos, 2*max_snp*sizeof(int));
      assert(pos != NULL);
      for (isamp = 0; isamp < nsample; isamp++) {
	geno[isamp] = realloc(geno[isamp], 2*max_snp*sizeof(int));
	assert(geno[isamp] != NULL);
      }
      freq = realloc(freq, 2*max_snp*sizeof(double*));
      assert(freq != NULL);
      for (isnp = max_snp; isnp < 2*max_snp; isnp++) {
	freq[isnp] = malloc((max_all+1) * sizeof(double));
	assert(freq[isnp] != NULL);
      }
      max_snp *= 2;
    }
    iline++;
    totall = killit = 0;
    for (iall = 0; iall <= max_all; iall++) {allcount[iall] = 0;}
    for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
      if (itoken == 0) {chr = strtol(token, NULL, 10);}
      else if (itoken == 1) {
	pos[nsnp] = strtol(token, NULL, 10);
	if ( (chr == prev_chrom && pos[nsnp] < pos[nsnp-1]) || chr < prev_chrom) {
	  fprintf(stderr, "Variants are out of order: chrom=%d previous chrom=%d pos=%d previous pos=%d\n",
		  chr, prev_chrom, pos[nsnp], pos[nsnp-1]);
	  exit(0);
	}
	if (nsnp > 0 && chr == prev_chrom && pos[nsnp] - pos[nsnp-1] < min_snp_sep) {
	  killit = 1;
	  break;
	}
      }
      else {
	all = strtol(token, NULL, 10);
	if (all > max_all) {
	  killit = 1;
	  ex_all++;
	  break;
	}
	geno[itoken-2][nsnp] = all;
	if (use_sample[itoken-2] == 1) {	  
	  if (all >= 0) {
	    allcount[all]++;
	    totall++;
	  }
	}
      }
    }
    // if reading freqs from file, read one
    if (freq_flag == 1) {
      fgets(newLine, linesize, ff);
      fpos = fchr = 0;
      for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
      	if (itoken == 0) {fchr = strtol(token, NULL, 10);}
      	else if (itoken == 1) {fpos = strtol(token, NULL, 10);}
      	else if (itoken > 1) {ffreq[itoken-2] = strtod(token, NULL);}
      }
      if (fchr != chr || fpos != pos[nsnp]) {
      	fprintf(stderr, 
	  "Mismatch between data file and frequency file. Data file (chr/pos): %d/%d vs freq file: %d/%d\n", 
      		chr, pos[nsnp], fchr, fpos);
      	abort();
      }
    }
    if (chr > nchrom) {continue;}

    // process this variant -- calculate allele frequencies
    maxall = 0;
    nall[nsnp] = 0;
    for (iall = 0; iall <= max_all; iall++) {
      if (allcount[iall] > 0) {
	nall[nsnp]++;
      }
      if (allcount[iall] > maxall) {
	maxall = allcount[iall];
	majall = iall;
      }
      if (freq_flag == 1) {
	freq[nsnp][iall] = ffreq[iall];
      }
      else {
	freq[nsnp][iall] = (double) allcount[iall] / totall;
      }
    }
    if (freq_flag == 1 && maxall == totall) {killit=1;}   // Skip monomorphic, unless external data on freq
    if (killit == 1) {continue;}
    prev_chrom = chr;
    if (nall[nsnp] > 2) {
      ntri++;
    }

    // Tabulate differences by pair for calculating discordance
    //   (I could skip this for use_pair == 0, but the saving is small)
    ipair = 0;
    for (isamp = 0; isamp < nsample; isamp++) {
      if (use_sample[isamp] == 0) {continue;}
      for (jsamp = isamp+1; jsamp < nsample; jsamp++) {
	if (use_sample[jsamp] == 0) {continue;}
	if (geno[isamp][nsnp] != -1 && geno[jsamp][nsnp] != -1) {
	  if (geno[isamp][nsnp] == geno[jsamp][nsnp]) {
	    if (geno[isamp][nsnp] != majall) {
	      same_min[ipair]++;
	    }
	  }
	  else {diff[ipair]++;}
	}
	else {
	  nmiss_bypair[ipair]++;
	}
	ipair++;
      }
    }

    if (nsnp < start_chr[chr]) {
      start_chr[chr] = nsnp;
    }
    if (nsnp > end_chr[chr]) {
      end_chr[chr] = nsnp;
    }
    nsnp++;
  }  // End loop over input snps

  // Weed out pairs with wrong discordance or too few informative markers
  fprintf(stdout, "%d variants,\t%d with >2 alleles\n", nsnp, ntri);
  ipair = 0;
  for (isamp = 0; isamp < nsample; isamp++) {
    if (use_sample[isamp] == 0) {continue;}
    for (jsamp = isamp+1; jsamp < nsample; jsamp++) {
      if (use_sample[jsamp] == 0) {continue;}
      sum = diff[ipair] + same_min[ipair];
      if (sum == 0) {
	discord[isamp][jsamp]  = 0;
      }
      else {
	discord[isamp][jsamp] = (double) diff[ipair] / sum;
      }
      if (use_pair[isamp][jsamp] == 1) {
	if (discord[isamp][jsamp] >  max_discord || sum < min_inform ||
	    discord[isamp][jsamp] < min_discord) {
	  use_pair[isamp][jsamp] = 0;
	}
	else {
	  nuse_pair++;
	}
      }
      ipair++;
    }
  }  
  fprintf(stdout, "sample pairs analyzed (filtered for discordance and informative markers): %d\n", nuse_pair);

  maxlen = 0;
  for (chr = 1; chr <= nchrom; chr++) {
    if (end_chr[chr] - start_chr[chr] + 1 > maxlen) {maxlen = end_chr[chr] - start_chr[chr] + 1;}
  }
  for (is = 0; is < 2; is++) {
    psi[is] = malloc(maxlen * sizeof(int));
    phi[is] = malloc(maxlen * sizeof(double));
    b[is] = malloc(maxlen * sizeof(double));
    alpha[is] = malloc(maxlen * sizeof(double));
    beta[is] = malloc(maxlen * sizeof(double));
  }
  scale = malloc(maxlen * sizeof(double));
  traj = malloc(maxlen * sizeof(int));

  ipair = 0;
  for (isamp = 0; isamp < nsample; isamp++) {
    if (use_sample[isamp] == 0) {continue;}
    for (jsamp = isamp+1; jsamp < nsample; jsamp++) {
      if (use_sample[jsamp] == 0) {continue;}
      sum = diff[ipair] + same_min[ipair];
      if (use_pair[isamp][jsamp] == 1) {
	last_prob = 0;
	last_pi = pi[0] = pinit[0];  // flat prior
	pi[1] = pinit[1];
	last_krec = k_rec = k_rec_init;
	finish_fit = 0;
	for (iter = 0; iter < niter; iter++) {
	  trans_obs = trans_pred = ntrans = 0;
	  seq_ibd = seq_dbd = count_ibd_fb = count_dbd_fb = seq_ibd_fb = seq_dbd_fb = 0;
	  max_phi = 0;
	  for (chr = 1; chr <= nchrom; chr++) {
	    if (end_chr[chr] == 0) {continue;}
	    chrlen = pos[end_chr[chr]] - pos[start_chr[chr]];
	    for (isnp = start_chr[chr]; isnp <= end_chr[chr]; isnp++) {
	      snp_ind = isnp - start_chr[chr];
	      gi = geno[isamp][isnp];
	      gj = geno[jsamp][isnp];
	      pright = 1 - eps * (nall[isnp] - 1);
	      if (gi == -1 || gj == -1) {
		// O = n (missing data)
		b[0][snp_ind] = b[1][snp_ind] = 0.5;
	      }
	      else if (gi == gj) {
		// homozygote
		b[0][snp_ind] = pright * pright * freq[isnp][gi] + eps * eps * (1 - freq[isnp][gi]); //IBD
		b[1][snp_ind] = pright * pright * freq[isnp][gi]*freq[isnp][gi] + 
		  2 * pright * eps * freq[isnp][gi] * (1-freq[isnp][gi]) + 
		  eps * eps * (1-freq[isnp][gi]) * (1-freq[isnp][gi]);
	      }
	      else {
		// O = aA
		b[0][snp_ind] = 2*pright*eps*(freq[isnp][gi]+freq[isnp][gj]) + 
		  2 * eps * eps * (1 - freq[isnp][gi] - freq[isnp][gj]);
		b[1][snp_ind] = 2*pright*pright*freq[isnp][gi]*freq[isnp][gj] + 
		  2*pright*eps*( freq[isnp][gi]*(1-freq[isnp][gj]) + freq[isnp][gj]*(1-freq[isnp][gi]) ) +
		  2 * eps * eps * (1 - freq[isnp][gi] - freq[isnp][gj] + freq[isnp][gi]*freq[isnp][gj]);
	      }
	      if (isnp == start_chr[chr]) {
		psi[0][snp_ind] = psi[1][snp_ind] = 0;
		for (is = 0; is < 2; is++) {
		  phi[is][snp_ind] = log(pi[is]) + log(b[is][snp_ind]);
		  alpha[is][snp_ind] = pi[is] * b[is][snp_ind];
		  scale[is] = 1.;
		}
	      }
	      else {
		ptrans = k_rec * rec_rate * (pos[isnp] - pos[isnp-1]);
		a[0][1] = 1 - pi[0] - (1 - pi[0]) * exp(-ptrans);
		a[1][0] = 1 - pi[1] - (1 - pi[1]) * exp(-ptrans);
		a[0][0] = 1 - a[0][1];
		a[1][1] = 1 - a[1][0];
		for (js = 0; js < 2; js++) {    // index over state of current snp
		  maxval = -10000000;
		  alpha[js][snp_ind] = scale[snp_ind] = 0;
		  for (is = 0; is < 2; is++) {    // index over state of previous snp
		    if (phi[is][snp_ind-1] + log(a[is][js]) > maxval ) {
		      maxval = phi[is][snp_ind-1] + log(a[is][js]);
		      psi[js][snp_ind] = is;
		    }
		    phi[js][snp_ind] = maxval + log(b[js][snp_ind]);
		    alpha[js][snp_ind] += alpha[is][snp_ind-1] * a[is][js] * b[js][snp_ind];
		  }
		  scale[snp_ind] += alpha[js][snp_ind];
		}
		for (js = 0; js < 2; js++) {    // scale alpha to prevent underflow (Rabiner eqns 92)
		  alpha[js][snp_ind] /= scale[snp_ind];
		}
	      }   // end if initializing/continuing
	    }  // end snp loop
	    last_snp = snp_ind;
	    max_phiL = phi[1][snp_ind];
	    if (phi[0][snp_ind] > phi[1][snp_ind]) {max_phiL = phi[0][snp_ind];}
	    max = (phi[1][snp_ind] > phi[0][snp_ind]) ? 1 : 0;
	    traj[snp_ind] = max;
	    max_phi += max_phiL;
	    
	    // Loop backward to calculate betas (backward part of forward-backward)
	    isnp = end_chr[chr];
	    snp_ind = isnp - start_chr[chr];
	    beta[0][snp_ind] = beta[1][snp_ind] = 1;
	    for (isnp = end_chr[chr]-1; isnp >= start_chr[chr]; isnp--) {
	      snp_ind = isnp - start_chr[chr];
	      ptrans = k_rec * rec_rate * (pos[isnp+1] - pos[isnp]);
	      a[0][1] = 1 - pi[0] - (1 - pi[0]) * exp(-ptrans);
	      a[1][0] = 1 - pi[1] - (1 - pi[1]) * exp(-ptrans);
	      a[0][0] = 1 - a[0][1];
	      a[1][1] = 1 - a[1][0];
	      for (is = 0; is < 2; is++) {    // index over state of current snp
		beta[is][snp_ind] = 0;
		for (js = 0; js < 2; js++) {    
		  // index over state of previous snp (= snp+1, since looping backward)
		  beta[is][snp_ind] += beta[js][snp_ind+1] * a[is][js] * b[js][snp_ind+1] / scale[snp_ind];
		}
	      }
	    }
	    // inverse loop for Viterbi
	    for (snp_ind = last_snp-1; snp_ind >= 0; snp_ind--) {
	      traj[snp_ind] = psi[max][snp_ind+1];
	      max = traj[snp_ind];
	    }
	    // tabulate FB results (could be done in previous loop, but I like clarity)
	    for (snp_ind = 0; snp_ind < last_snp; snp_ind++) {
	      p_ibd = alpha[0][snp_ind] * beta[0][snp_ind] / 
		(alpha[0][snp_ind] * beta[0][snp_ind] + alpha[1][snp_ind] * beta[1][snp_ind]);
	      count_ibd_fb += p_ibd;
	      count_dbd_fb += (1-p_ibd);
	      if (snp_ind < last_snp-1) {
		isnp = snp_ind + start_chr[chr];
		delpos = pos[isnp+1] - pos[isnp];
		ptrans = k_rec * rec_rate * delpos;
		a[0][1] = 1 - pi[0] - (1 - pi[0]) * exp(-ptrans);
		a[1][0] = 1 - pi[1] - (1 - pi[1]) * exp(-ptrans);
		a[0][0] = 1 - a[0][1];
		a[1][1] = 1 - a[1][0];
		xisum = 0;
		for (is = 0; is < 2; is++) {
		  for (js = 0; js < 2; js++) {
		    xi[is][js] = alpha[is][snp_ind] * a[is][js] * b[js][snp_ind+1] * beta[js][snp_ind+1];
		    xisum += xi[is][js];
		  }
		}
		for (is = 0; is < 2; is++) {
		  gamma[is] = 0;
		  for (js = 0; js < 2; js++) {
		    xi[is][js] /= xisum;
		    gamma[is] += xi[is][js];
		  }
		}
		// xi(0,0) in Rabiner notation = prob of being in state 0 at t and 0 at t+1
		seq_ibd_fb += delpos * xi[0][0];
		seq_dbd_fb += delpos * xi[1][1];
		trans_obs += (xi[0][1] + xi[1][0]);
		trans_pred += gamma[0] * a[0][1] + gamma[1] * a[1][0];
	      }
	    }
	    
	    // print final Viterbi trajectory
	    // start
	    if (iter == niter - 1 || finish_fit != 0) {
	      fprintf(outf, "%s\t%s\t%d\t%d", sample[isamp], sample[jsamp], chr, pos[0+start_chr[chr]]);
	    }
	    start_pos = pos[0+start_chr[chr]];
	    for (isnp = 1; isnp < end_chr[chr] - start_chr[chr] + 1; isnp++) {
	      if (traj[isnp] != traj[isnp-1]) {
		ntrans++;
		add_seq = pos[isnp - 1 + start_chr[chr]] - start_pos + 1;
		if (traj[isnp-1] == 0) {seq_ibd += add_seq;}
		else {seq_dbd += add_seq;}
		start_pos = pos[isnp + start_chr[chr]];
		// end one and start another
		if (iter == niter -1 || finish_fit != 0) {
		  fprintf(outf, "\t%d\t%d\n", pos[isnp - 1 + start_chr[chr]], traj[isnp-1]);
		  fprintf(outf, "%s\t%s\t%d\t%d", sample[isamp], 
			  sample[jsamp], chr, pos[isnp + start_chr[chr]]);
		}
	      }
	    }
	    isnp = end_chr[chr] - start_chr[chr];
	    add_seq = pos[isnp + start_chr[chr]] - start_pos + 1;
	    if (traj[isnp] == 0) {seq_ibd += add_seq;}
	    else {seq_dbd += add_seq;}
	    if (iter == niter -1 || finish_fit != 0) {
	      fprintf(outf, "\t%d\t%d\n", pos[end_chr[chr]], traj[isnp]);
	    }
	  }  // end chrom loop
	  
	  // quit if fit converged on previous iteration; otherwise, update parameters
	  if (finish_fit != 0) {break;}
	  pi[0] = count_ibd_fb / (count_ibd_fb + count_dbd_fb);
	  k_rec *= trans_obs / trans_pred;
	  // ad hoc attempt to avoid being trapped in extremum
	  if (iter < niter-1 && finish_fit == 0) {
	    if (pi[0] < 1e-5) {pi[0] = 1e-5;}
	    else if (pi[0] > 1- 1e-5) {pi[0] = 1 - 1e-5;}
	    if (k_rec < 1e-5) {k_rec = 1e-5;}
	  }
	  pi[1] = 1 - pi[0];
	  delpi = pi[0] - last_pi;
	  delk = k_rec - last_krec;
	  delprob = max_phi - last_prob;
	  last_pi = pi[0];
	  last_krec = k_rec;
	  last_prob = max_phi;

	  // Evaluate fit
	  if (fabs(delpi) < fit_thresh_dpi && 
	       (fabs(delk) < fit_thresh_dk || fabs(delk/k_rec) < fit_thresh_drelk) ) {
	    finish_fit = 1;
	  }
	}  // end parameter fitting loop

	// probabilistic shared seq output: 
	// seq_ibd_fb / (seq_ibd_fb + seq_dbd_fb), 
	fprintf(pf, "%s\t%s\t%d\t%.4f\t%.5e\t%d\t%.3f\t%d\t%.5f\t%.5f\n", 
		sample[isamp], sample[jsamp], sum, discord[isamp][jsamp],
		max_phi, iter, k_rec, ntrans, (double) seq_ibd / (seq_ibd + seq_dbd),
		count_ibd_fb / (count_ibd_fb + count_dbd_fb));
      }   // end if use pair
      ipair++;
      if (ipair%1000 == 0) {fprintf(stdout, "Starting pair %d\n", ipair);}
    }
  }
  if (ex_all > 0) {
    fprintf(stdout, "Variants with too many alleles: %d\n", ex_all);
  }
  return(0);
}
