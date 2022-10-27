import argparse
import numpy as np
from pomegranate import *
import os

def hidden_states(bsTot_list,bsMeth_list):
    ratio = np.divide(bsMeth_list,bsTot_list)
    case = np.nanmean(ratio[:,l1_list],axis=1)
    control = np.nanmean(ratio[:,l2_list],axis=1)
    diff = case - control
    return diff

def get_counts():
        tot_list = {}
        meth_list = {}
        for bsTot_file,bsMeth_file in zip(total_files,methylated_files):
                bsTot_list = np.loadtxt(bsTot_file,delimiter='\t',skiprows=0,dtype='int', usecols=range(1, num_samples+1))
                bsMeth_list = np.loadtxt(bsMeth_file,delimiter='\t',skiprows=0,dtype='int', usecols=range(1, num_samples+1))
                loci = np.genfromtxt(bsTot_file,delimiter='\t',skip_header=0, dtype='str', usecols=(0,))
                for x,y,z in zip(loci,bsTot_list,bsMeth_list):
                        tot_list[x] = y
                        meth_list[x] = z
        return(tot_list,meth_list)

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='LuxPom')
        parser.add_argument('-l1', '--l1', action='store', dest='l1', type=str, required=True, help='comma-delimited list containing control indices, starting from zero')
        parser.add_argument('-l2', '--l2', action='store', dest='l2', type=str, required=True, help='comma-delimited list containing case indices, starting from zero')
        parser.add_argument('-d1', '--data_total', action='store', dest='data_total', type=str, required=True, help='file containing input textfiles for total read counts; each line contains the file which contains data for one chromosome')
        parser.add_argument('-d2', '--data_methylated', action='store', dest='data_methylated', type=str, required=True, help='file containing input textfiles for methylated read counts; each line contains the file which contains data for one chromosome')
        parser.add_argument('-o','--outfolder', action='store', dest='outfolder', type=str, required=False, default='hmm_output', help='output location')
        parser.add_argument('-m','--min_count', action='store', dest='min_count', type=int, required=False, default=5, help='minimum total read count')
        parser.add_argument('-c','--min_CpGs', action='store', dest='min_CpGs', type=int, required=False, default=2, help='minimum number of CpGs in regions')
        parser.add_argument('-v','--version',action='version',version='0.666')
        options = parser.parse_args()

        l1_list = list(map(int,options.l1.split(',')))
        l2_list = list(map(int,options.l2.split(',')))
        total_files = [infile.strip() for infile in open(options.data_total,'r').readlines()]
        methylated_files = [infile.strip() for infile in open(options.data_methylated,'r').readlines()]
        outfolder = options.outfolder
        if not os.path.exists(outfolder): os.makedirs(outfolder)
    
        num_samples = len(l1_list) + len(l2_list)
        diffs = []
        for bsTot_file,bsMeth_file in zip(total_files,methylated_files):
                bsTot_list = np.loadtxt(bsTot_file,delimiter='\t',skiprows=0,dtype='float', usecols=range(1, num_samples+1))
                bsTot_list[bsTot_list < options.min_count] = np.nan
                bsMeth_list = np.loadtxt(bsMeth_file,delimiter='\t',skiprows=0,dtype='float', usecols=range(1, num_samples+1))
                diff=hidden_states(bsTot_list,bsMeth_list)
                diffs.append(diff)
            
        s0 = State(NormalDistribution(0, 0.08), name = 's0')
        s1 = State(NormalDistribution(0.3, 0.06), name = 's1')
        s2 = State(NormalDistribution(-0.3, 0.06), name = 's2')
        model = HiddenMarkovModel()
        model.add_states(s0, s1, s2)
        model.add_transition(model.start, s0, 0.333)
        model.add_transition(model.start, s1, 0.333)
        model.add_transition(model.start, s2, 0.333)
        model.add_transition(s0, s0, 0.5)
        model.add_transition(s0, s1, 0.25)
        model.add_transition(s0, s2, 0.25)
        model.add_transition(s1, s0, 0.25)
        model.add_transition(s1, s1, 0.5)
        model.add_transition(s1, s2, 0.25)
        model.add_transition(s2, s0, 0.25)
        model.add_transition(s2, s1, 0.25)
        model.add_transition(s2, s2, 0.5)
        model.add_transition(s0, model.end, 0.333)
        model.add_transition(s1, model.end, 0.333)
        model.add_transition(s2, model.end, 0.333)
        model.bake()
        model.fit(diffs, distribution_inertia=1.0)
        outfile = '%s/model_states.txt'%outfolder
        fo = open(outfile,'w')
        fo.write('%s'%model.states)
        fo.close()
        outfile = '%s/transition_probs.txt'%outfolder
        fo = open(outfile,'w')
        fo.write('%s'%model.dense_transition_matrix())
        fo.close()

        outfile = '%s/hidden_states_all.txt'%(outfolder)
        fo = open(outfile,'w')
        for i in range(len(total_files)):
	        bsTot_file = total_files[i]
	        loci = np.genfromtxt(bsTot_file,delimiter='\t',skip_header=0, dtype='str', usecols=(0,))
	        diff = diffs[i]
	        viterbi_likelihood, viterbi_path = model.viterbi(diff)
	        for r,s,t in zip(diff,viterbi_path[1:],loci):
		        fo.write('%s\t%s\t%s\n'%(t,s[1].name,r))
        fo.close()


        tot_list,meth_list = get_counts()
        infile = '%s/hidden_states_all.txt'%(outfolder)
        lines = open(infile,'r').readlines()
        outfile1 = '%s/total_reads_all.txt'%(outfolder)
        fo1 = open(outfile1,'w')
        outfile2 = '%s/methylated_reads_all.txt'%(outfolder)
        fo2 = open(outfile2,'w')
        outfile3 = '%s/counts.txt'%(outfolder)
        fo3 = open(outfile3,'w')

        lines = open(infile,'r').readlines()
        line = lines[0].strip().split()
        prev_state = line[1]
        prev_tot_count = tot_list[line[0]]
        prev_meth_count = meth_list[line[0]]
        count = 1
        prev_chrom, coord = line[0].split(':')
        start = coord
        end = coord
        for line in lines[1:]:
                line = line.strip().split()
                chrom, coord = line[0].split(':')
                tot_count = tot_list[line[0]]
                meth_count = meth_list[line[0]]
                state = line[1]
                if chrom == prev_chrom:
                        if state == prev_state:
                                count += 1
                                prev_tot_count += tot_count
                                prev_meth_count += meth_count

                        else:
                                if count > options.min_CpGs and prev_state != 's0':
                                        fo1.write('%s:%s:%s\t%s\n'%(chrom,start,end,'\t'.join(prev_tot_count.astype(str))))
                                        fo2.write('%s:%s:%s\t%s\n'%(chrom,start,end,'\t'.join(prev_meth_count.astype(str))))
                                fo3.write('%s:%s:%s\t%s\t%s\n'%(chrom,start,end,prev_state,count))
                                count = 1
                                start = coord
                                prev_tot_count, prev_meth_count = tot_count, meth_count
                else:
                        if count > options.min_CpGs and prev_state != 's0':
                                fo1.write('%s:%s:%s\t%s\n'%(chrom,start,end,'\t'.join(prev_tot_count.astype(str))))
                                fo2.write('%s:%s:%s\t%s\n'%(chrom,start,end,'\t'.join(prev_meth_count.astype(str))))
                        fo3.write('%s:%s:%s\t%s\t%s\n'%(chrom,start,end,prev_state,count))
                        count = 1
                        start = coord
                        prev_tot_count, prev_meth_count = tot_count, meth_count
                prev_state = state
                prev_chrom = chrom
                end = coord
        if count > options.min_CpGs and prev_state != 's0':
                fo1.write('%s:%s:%s\t%s\n'%(chrom,start,end,'\t'.join(prev_tot_count.astype(str))))
                fo2.write('%s:%s:%s\t%s\n'%(chrom,start,end,'\t'.join(prev_meth_count.astype(str))))
        fo3.write('%s:%s:%s\t%s\t%s\n'%(chrom,start,end,prev_state,count))
        fo1.close()
        fo2.close()
        fo3.close()

