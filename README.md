# LuxPom

## Overview
LuxRep is a genome wide methylation analysis tool that models different bisulphite-converted DNA libraries from the same biological replicate (i.e., technical replicates) thereby allowing use of libraries of even low bisulphite conversion rates. Faster processing suited for genome wide analysis is achieved by using variational inference for posterior approximation. These features extend on **LuxGLM**, a probabilistic method for methylation analysis that handles complex sample designs.

## Features

* Model-based integration of biological and technical replicates
* Full Bayesian inference by variational inference implemented in __Stan__

## Quick introduction

An usual LuxRep pipeline has the following steps

1. Generate count data from sequencing reads using e.g. **Bismark**
	1. Align BS-seq data
	2. Extract converted and unconverted counts
2. Format input files (sample scripts using **bedtools** provided below)
3. Methylation analysis with **LuxRep** (instructions given below)
	1. Estimate experimental parameters from control cytosine data
	2. Quantify methylation modification and identify differential methylation

## Installation

Below the environmental parameter $STAN_HOME refers to the source directory of **CmdStan**

LuxRep requires the following software

* __CmdStan__ (tested on version 2.18.0)
* __Python__ (tested on version 3.7.5)
* __pystan__ (tested on version 2.17.1.0)
* __Numpy__ (tested on version 1.20.2)
* __Scipy__ (tested on version 1.1.0)

In addition, data preprocessing can be facilitated by the following

* __bedtools__ (tested on version 2.27.1)

## Using LuxPom

The first step in LuxPom is to segment the genome via HMM (using the python package pomegranate). A python script **run\_hmm.py** is supplied which calls the python package `pomegranate` to determine the methylation state for each cytosine then combines adjacent cytosines with the same methylation state into regions (includes compiling of relevant __Stan__ code). The script returns the output files **total_reads_all.txt** and **methylated_reads_all.txt** which contain the total and methylated reads, respectively, for the segmented regions.

	 usage: run_hmm.py [-h] -l1 CONTROL_INDICES_LIST -l2 CASE_INDICES_LIST -d1 TOTAL_READ_COUNTS_FILE -d2 METHYLATED_READ_COUNTS_FILE -o OUTFOLDER -m MIN_TOTAL_COUNT -c MIN_CPGS

	 
	 Estimates experimental parameters bsEff and seqErr
	 
	 optional arguments:
	 -h, --help									show this help message and exit
	 -l1 CONTROL_INDICES_LIST, --l1 CONTROL_INDICES_LIST				comma-delimited list containing control indices, starting from zero
	 -l2 CASE_INDICES_LIST, --l2 CASE_INDICIES_LIST					comma-delimited list containing case indices, starting from zero
	 -d1 TOTAL_READ_COUNTS_FILE, --data_total TOTAL_READ_COUNTS_FILE		file containing input textfiles for total read counts; one line per textfile
	 -d2 METHYLATED_READ_COUNTS_FILE, --data_methylated METHYLATED_READ_COUNTS_FILE	file containing input textfiles for methylated read counts; one line per textfile
	 -o OUTFOLDER, --outfolder OUTFOLDER 						output location
	 -m MIN_TOTAL_COUNT, --min_count MIN_TOTAL_COUNT				minimum total read count; if not specified, 5 is used as default
	 -c MIN_CPGS, --min_CpGs MIN_CPGS						minimum number of CpGs in regions; if not specified, 2 is used as default
	 -v, --version									show program's version number and exit

For instance, run\_hmm.py can be called as
	
	python run_hmm.py run_hmm.py -l1 0,1,2,3,4,5 -l2 6,7,8,9,10,11 -c 5 -d1 total_files.txt -d2 methylated_files.txt
	python run_hmm.py -f data/lambda_fileList.txt -n data/nameList.txt -o $PWD/controls -l $STAN_HOME

*Input*

The file **lambda\_fileList.txt** contains a space-delimited list of bed files containing read counts of methylated and unmethylated cytosines generated using the **Bismark** pipeline. Each tab delimited bed file follows the format `<chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>`. The sample bed files in this repo contain '-' in the fourth column in place of methylation percentage.

The file **nameList.txt** contains a space-delimited list of labels for the libraries corresponding to the files in **lambda\_fileList.txt**.

*Output*

The output files consist of __Stan__-generated diagnostic and output files including **output.csv** which contains samples from the approximate posterior of the model parameters. The mean posterior for the parameters bsEff and seqErr from the latter are used as estimates of the experimental parameters used in the second module. 

The second step in using LuxRep is estimating the methylation levels of the noncontrol cytosine and testing for differential methylation. A python script **luxrep.py** is supplied for generating input files from user-supplied data files and running the analysis (includes compiling of relevant __Stan__ code). 

	 usage: luxrep.py -d NONCONTROL_DATA -s SAMPLE_LIST -m DESIGN_MATRIX -n LIBRARY_LABELS_LIST -c EXPERIMENTAL_PARAMETERS -l $STAN_HOME
	 
	 Estimates experimental parameters bsEff and seqErr
	 
	 optional arguments:
	 -h, --help								show this help message and exit
	 -d NONCONTROL_DATA, --data NONCONTROL_DATA				file containing noncontrol cytosine data
	 -s SAMPLE_LIST, --sample_list SAMPLE_LIST				file containing sample number of libraries
	 -m DESIGN_MATRIX, --design_matrix					file containing design matrix
	 -n LIBRARY_LABELS_LIST, --library_names LIBRARY_LABELS_LIST		space-delimited list of labels for files listed in lambda_fileList.txt
	 -c EXPERIMENTAL_PARAMETERS, --exp_params EXPERIMENTAL_PARAMETERS	directory containing output for control data with full pathname
	 -o OUTFOLDER, --outfolder OUTFOLDER					directory containing data analysis output with full pathname
	 -l $STAN_HOME, --cmdstan_loc $STAN_HOME 				CmdStan directory with full pathname
	 -v, --version								show program's version number and exit

For instance, luxrep.py can be called as

	python luxrep.py -d results/counts_1.tab -s data/sample_list.txt -m data/design_matrix.txt -n data/nameList.txt -c $PWD/controls -o $PWD/results/1 -l $STAN_HOME

*Input*

The file **counts_1.tab** contains noncontrol cytosine data for all libraries in the dataset with data for each library corresponding to a two-column block, (i) total read count and (ii) "C" count. The first row contains labels for the libraries and each row thereafter contains data for one cytosine (this table format is the same as the one used by the methylation analysis tool **Radmeth**). For M libraries, each row contains 1 + 2xM tab-delimited columns with the first column showing the genomic coordinates (in the format `<chromosome>:<start>:<end>`) and the rest showing N<sub><sub>BS</sub>1</sub>, N<sub><sub>BS,C</sub>1</sub>, N<sub><sub>BS</sub>2</sub>, N<sub><sub>BS,C</sub>2</sub>, ..., N<sub><sub>BS</sub>M</sub>, N<sub><sub>BS,C</sub>M</sub>. The order of libraries should correspond with the order in the data files for the control cytosine used in the previous module. 

When using output coverage files from **Bismark**'s pipeline (see file format description above in section for control data), with the **bedtools** toolset the following bash scripts may be used to merge the coverage files into the format required by **luxrep.py**. The file **fileList.txt** is a space-delimited list of the bed files (one bed file per library) and **nameList.txt** is a space-delimited list of library labels.

	# cd into the folder containing the bed files, generate intermediate files, one for each of the M libraries
	for i in {1..M}; do file=$i".bed"; ls $file; rm -f $i"_temp.bed"; cat $file | awk '{print $1 "\t" $2 "\t" $3 "\t" $5+$6 ";" $5}' > $i"_temp.bed"; done
	
	# combine the intermediate bed files
	unionBedGraphs -i $(cat ../fileList.txt) -names $(cat ../nameList.txt) -header -filler "0;0" | sed 's/'$'\t/:/;s/'$'\t/:/;s/;/'$'\t/g;s/chrom:start:end'$'\t//g' > ../hg19.tab

The code snippets above loops over the bed files and for each, sums columns 5 and 6 (methylated and unmethylated read counts, respectively) and combines it with column 5 into one column (separated by a semicolon) to comply with the format requirement of **unionBedGraphs**. The output file from **unionBedGraphs** is further formatted to combine the genomic coordinates into one column and expand the total and methylated counts into two columns, etc.

*Optional preprocessing*

Taken from **Radmeth**'s pipeline, an optional step allows pseudo-parallelization of running the methylation analysis module by splitting the input table into smaller tables that can be analysed independently. A python script is provided for this step.

	 usage: split_inputTable.py -d NONCONTROL_DATA -s SPLIT_SIZE -o OUTFOLDER
	 
	 Splits coverage table
	 
	 optional arguments:
	 -h, --help					show this help message and exit
	 -d NONCONTROL_DATA, --data NONCONTROL_DATA	file containing noncontrol cytosine data
	 -s SPLIT_SIZE, --split_size SPLIT_SIZE		number of cytosine in split table
	 -o OUTFOLDER, --outfolder OUTFOLDER		directory containing data analysis output with full pathname
	 -v, --version					show program's version number and exit

For instance, split\_inputTable.py can be called as

	python split_inputTable.py -d data/hg19.tab -s 1000 -o $PWD/results

which generates several files each containing data for <= 1000 cytosine.

The file **sample\_list.txt** lists the sample membership ("biological replicate") of each of the M libraries ("technical replicates"). The file contains one column and M rows in the same order of the libraries in the data files for the control and noncontrol cytosine (i.e., library<sub>row1</sub> = library<sub>column1</sub>). Each row contains one integer x<sub>_i_</sub> &in; 1 ... N<sub>sample</sub>, where _i_ &in; 1 ... M.

The file **design\_matrix.txt** contains covariate information (like in **LuxGLM**) but with one additional row on top and one column on the left for covariate and sample labels (biological replicates), respectively. The row order should follow the sample order in **sample\_list.txt** with N<sub>rows</sub> = N<sub>sample</sub>. The default order of the columns, after the row header, starts with the intercept ("base") followed by the covariate to be tested for differential methylation (in the given example, the covariate "case"), then, if present, followed by the other covariates.  

The directory specified by EXPERIMENTAL\_PARAMETERS contains the output files from the first module. The second module looks for the file **output.csv** and parses it for the posterior mean of the parameters bsEff and seqErr. 

*Output*

The directory specified by OUTFOLDER contains N<sub>cytosine</sub> subfolders each containing input, diagnostic and output files with the latter showing posterior samples of the model parameters for the corresponding cytosine. The second module includes a routine for computing a Bayes factor for differential methylation based on the covariate specified in **design\_matrix.txt**. A summary text file (in **.bed** format) lists the genomic coordinates for each cytosine in the input data file followed by its Bayes factor (**bfs.bed**).

To test for significance of a covariate other than the default (second column in the design matrix), a routine **savagedickey2** in the library file **luxrep_routines.py** is provided. The routine requires a text file containing a list of **output.csv** files with full pathname (one per line) for all cytosine to be analysed and a number corresponding to the column in the design matrix of the covariate of interest (see **TUTORIAL**). The routine generates a bed file in the current directory with a list of Bayes factor for each cytosine in the same order as the input text file.

**References**

[1] T. Äijö, X. Yue, A. Rao and H. Lähdesmäki, “LuxGLM: a probabilistic covariate model for quantification of DNA methylation modifications with complex experimental designs.,” Bioinformatics, 32.17:i511-i519, Sep 2016.

[2] F. Krueger and S. R. Andrews, “Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications.,” Bioinformatics, 27.11:1571-1572, Jun 2011. 

[3] Stan Development Team, “Stan Modeling Language Users Guide and Reference Manual, Version 2.14.0.,” http://mc-stan.org, Jan 2016. 

[4] Stan Development Team, “PyStan: the Python interface to Stan, Version 2.14.0.0.,” http://mc-stan.org, Jan 2016.

[5] A. R. Quinlan and I. M. Hall, “BEDTools: a flexible suite of utilities for comparing genomic features,” Bioinformatics, 26.6:841–842, Mar 2010.

