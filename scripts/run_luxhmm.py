#!/usr/bin/env python

import argparse
import cmdstanpy
import sys
import numpy
import numpy.random
import scipy.stats
import os
import shutil

import random
from glob import glob
import stan
import pickle
import logging
from hashlib import md5

import luxhmm_routines

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='LuxHMM')
        parser.add_argument('-t', '--total_file', action='store', dest='total_file', type=str, required=False, default='hmm_output/total_reads_all.txt', help='file containing total read counts in regions; if not supplied, defaults to hmm_output/total_reads_all.txt')
        parser.add_argument('-m', '--methylated_file', action='store', dest='methylated_file', type=str, required=False, default='hmm_output/methylated_reads_all.txt', help='file containing methylated read counts in regions; if not supplied, defaults to hmm_output/methylated_reads_all.txt')
        parser.add_argument('-d', '--design_matrix', action='store', dest='design_file', type=str, required=True, help='file containing design matrix')
        parser.add_argument('-o', '--outfolder', action='store', dest='outfolder', type=str, required=False, default='{}/results'.format(os.getcwd()), help='directory containing data analysis output; if not supplied, defaults to results')
        parser.add_argument('-a', '--bsEff', action='store', dest='bsEff', type=str, required=False, help='comma-delimited bisulfite conversion efficiencies of replicates; defaults to 1.0')
        parser.add_argument('-b', '--bsBEff', action='store', dest='bsBEff', type=str, required=False, help='comma-delimited incorrect bisulfite conversion efficiencies of replicates; defaults to 0')
        parser.add_argument('-c', '--seqErr', action='store', dest='seqErr', type=str, required=False, help='comma-delimited sequencing error rates of replicates; defaults to 0')
        parser.add_argument('-r', '--run_advi', action='store', dest='advi', type=str, required=False, default='T', help='whether to run ADVI or not (HMC); if not supplied, defaults to T')
        parser.add_argument('-l','--cmdstan_loc', action='store', dest='cmdstan_directory', type=str, required=False, default='/scratch/cs/csb/users/malonzm1/software/cmdstan-2.29.0', help='cmdstan directory with full pathname')
        parser.add_argument('-v','--version',action='version',version='0.666')
        options = parser.parse_args()

        with open(options.design_file) as f:
                ncols = len(f.readline().split())

        D = numpy.loadtxt(options.design_file,delimiter='\t',skiprows=1,dtype='float', usecols=range(1, ncols + 1))
        nrows = D.shape[0]
        ncols2 = nrows
        # read data files
        if  options.bsEff is None:
                bsEff=numpy.ones(nrows)*1

        else:
                bsEff=list(map(float,options.bsEff.split(',')))

        if  options.bsBEff is None:
                bsBEff=numpy.ones(nrows)*0

        else:
                bsBEff=list(map(float,options.bsBEff.split(',')))

        if  options.seqErr is None:
                seqErr=numpy.ones(nrows)*0

        else:
                seqErr=list(map(float,options.seqErr.split(',')))

        bsTot_file=options.total_file
        bsTot_list = numpy.loadtxt(bsTot_file,delimiter='\t',skiprows=0,dtype='int', usecols=range(1, ncols2 + 1))
        if len(bsTot_list.shape) == 1:
                bsTot_list = numpy.array([bsTot_list])
        bsMeth_file=options.methylated_file
        bsMeth_list = numpy.loadtxt(bsMeth_file,delimiter='\t',skiprows=0,dtype='int', usecols=range(1, ncols2 + 1))
        if len(bsMeth_list.shape) == 1:
                bsMeth_list = numpy.array([bsMeth_list])
        loci = numpy.genfromtxt(bsTot_file,delimiter='\t',skip_header=0, dtype='str', usecols=(0,))
        if len(loci.shape) == 0:
                loci = numpy.array([loci])

        outfolder = options.outfolder
        if not os.path.exists(outfolder): os.makedirs(outfolder)

        stan_dir = options.cmdstan_directory

        curdir = os.getcwd()
        outfile1 = '{}/regions.bed'.format(outfolder)
        os.system('rm -f {}'.format(outfile1))
        os.system('touch {}'.format(outfile1))
        os.chdir(stan_dir)
        os.system('make {}/luxhmm'.format(curdir))
        for n in range(len(loci)):
                bsTot = bsTot_list[n]
                bsMeth = bsMeth_list[n]
                locus = loci[n]
            # get data and init dictionaries for stan
                data = luxhmm_routines.get_stan_input(bsTot,bsMeth,D,bsEff,bsBEff,seqErr)
                print(n+1)
                outdir = '{}/{}'.format(outfolder, locus)
                if not os.path.exists(outdir): os.makedirs(outdir)
                cmdstanpy.write_stan_json('{}/data.json'.format(outdir), data)

                os.chdir(outdir)
                os.system('cp {}/luxhmm .'.format(curdir))
                os.system('chmod 755 luxhmm')
                if options.advi == 'F': os.system('./luxhmm sample num_samples=250 num_chains=4 data file=data.json output file=output.csv')
                if options.advi == 'T': os.system('./luxhmm variational output_samples=1000 elbo_samples=1000 grad_samples=10 data file=data.json output diagnostic_file=diagnostics.csv > summary.txt')
                if options.advi == 'T': 
                        while 'COMPLETED' not in open('summary.txt','r').readlines()[-1]: os.system('./luxhmm variational output_samples=1000 elbo_samples=1000 grad_samples=10 data file=data.json output diagnostic_file=diagnostics.csv > summary.txt')
                # compute bayes factor  

                if options.advi == 'T': bf=luxhmm_routines.savagedickey_advi(locus,ncols)
                if options.advi == 'F': bf=luxhmm_routines.savagedickey_hmc(locus,ncols)
                #bf=luxhmm_routines.savagedickey(locus,ncols)
                fo1 = open(outfile1,'a')
                fo1.write('{}\t{}\n'.format('\t'.join(locus.split(':')),bf))
                fo1.close()
                shutil.rmtree(outdir)

