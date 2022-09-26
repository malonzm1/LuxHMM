# LuxPom

## Overview
LuxPom is a genome wide methylation analysis tool that detects differentially methylated regions from bisulfite sequencing data. It performs genome segmentation of cytosines into methylated regions using hidden Markov model (HMM) as implemented in the Python package `pomegranate` and infers differential methylation by logistic regression based on **LuxGLM** and **LuxUS**, probabilistic methods for methylation analysis that handle complex experimental deisgns. The model is implemented in __Stan__ and faster processing suited for genome wide analysis is achieved by using variational inference for posterior approximation, as featured in __Stan__. 

## Features

* Genome segmentation into methylated regions by HMM as implemented in the Python package `pomegranate`
* Full Bayesian inference by variational inference implemented in __Stan__

## Quick introduction

A usual LuxPom pipeline has the following steps

1. Generate count data from sequencing reads using e.g. **Bismark**
	1. Align BS-seq data
	2. Extract converted and unconverted counts
2. Methylation analysis with **LuxPom** (instructions given below)
	1. Segment genome into methylated regions
	2. Quantify methylation levels and identify differentially methylation regions

## Installation

Below the environmental parameter $STAN_HOME refers to the source directory of **CmdStan**

LuxRep requires the following software

* __CmdStan__ (tested on version 2.29.0)
* __pomegranate__ (tested on version 0.14.8)
*  __cmdstanpy__ (tested on version 1.0.1)
* __Python__ (tested on version 3.8.12)
* __pystan__ (tested on version 3.1.1)
* __Numpy__ (tested on version 1.22.2)
* __Scipy__ (tested on version 1.8.0)


## Using LuxPom

The first step in LuxPom is to segment the genome via HMM (using the python package `pomegranate`) into regions that are hypomethylated, hypermethylated or equally methylated between two groups. A python script **run\_hmm.py** is supplied which calls the python package `pomegranate` to determine the methylation state for each cytosine then combines adjacent cytosines with the same methylation state into regions. The script returns the output files **total_reads_all.txt** and **methylated_reads_all.txt** which contain the total and methylated reads, respectively, for the segmented regions.

	 usage: run_hmm.py [-h] -l1 CONTROL_INDICES_LIST -l2 CASE_INDICES_LIST -d1 TOTAL_READ_COUNTS_FILE -d2 METHYLATED_READ_COUNTS_FILE -o OUTFOLDER -m MIN_TOTAL_COUNT -c MIN_CPGS

	 
	 Estimates experimental parameters bsEff and seqErr
	 
	 optional arguments:
	 -h, --help									show this help message and exit
	 -l1 CONTROL_INDICES_LIST, --l1 CONTROL_INDICES_LIST				comma-delimited list containing control indices, starting from zero
	 -l2 CASE_INDICES_LIST, --l2 CASE_INDICIES_LIST					comma-delimited list containing case indices, starting from zero
	 -d1 TOTAL_READ_COUNTS_FILE, --data_total TOTAL_READ_COUNTS_FILE		file containing input textfiles for total read counts; each line contains the file which contains data for one chromosome
	 -d2 METHYLATED_READ_COUNTS_FILE, --data_methylated METHYLATED_READ_COUNTS_FILE	file containing input textfiles for methylated read counts; each line contains the file which contains data for one chromosome
	 -o OUTFOLDER, --outfolder OUTFOLDER 						output location; if not supplied, defaults to hmm_output
	 -m MIN_TOTAL_COUNT, --min_count MIN_TOTAL_COUNT				minimum total read count; if not specified, 5 is used as default
	 -c MIN_CPGS, --min_CpGs MIN_CPGS						minimum number of CpGs in regions; if not specified, 2 is used as default
	 -v, --version									show program's version number and exit

For instance, run\_hmm.py can be called as
	
	python run_hmm.py run_hmm.py -l1 0,1,2,3,4,5 -l2 6,7,8,9,10,11 -c 5 -d1 total_fileList.txt -d2 methylated_fileList.txt

*Input*

The file **total\_fileList.txt** contains the list of textfiles containing the total read counts, one file per line. Each textfile in the list contains the total read counts for one chromosome. The order of the columns must be consistent with that of CONTROL_INDICES_LIST and CASE_INDICES_LIST. Each line in the textfile is a tab-delimited list of total read counts, one line per CpG (must be sorted by chromosomal position). Each tab-delimited file follows the format `<chromosome>:<CpG position> <Total count, sample 1> ... <Total count, sample N>`. The file **methylated\_fileList.txt** contains the list of textfiles containing the methylated read counts, one file per line. Each textfile in the list contains the methylated read counts for one chromosome. Each line in the textfile is a tab-delimited list of methylated read counts, one line per CpG. Each tab-delimited file follows the format `<chromosome>:<CpG position> <Methylation count, sample 1> ... <Methylation count, sample N>`. 

*Output*

The output files are **total\_reads\_all.txt**, **methylated\_reads\_all.txt**, **hidden\_states\_all.txt**, **counts.txt**, **model\_states.txt** and **transition\_probs.txt**.

* **total\_reads\_all.txt** and **methylated\_reads\_all.txt** contain the combined total and methylated read counts, respectively, with one line per region. Each tab-delimited file follows the format `<chromosome>:<start>:<end> <Total/methylated count, sample 1> ... <Total/methylated count, sample N>`. These files are used as input for **run\_luxPom.py** described below.
* **hidden\_states\_all.txt** contains the methylation state for each CpG (s0 - equal methylation, s1 - hypermethylation and s2 - hypomethylation).
* **counts.txt** contains the methylation state and the number of CpGs for each region.
* **model\_states.txt** contains the state distributions.
* **transition\_probs.txt** contains the transition probabilities.


The second step in using LuxPom is estimating the methylation levels and inferring differential methylation. A python script **run\_luxPom.py** is supplied for generating input files from user-supplied data files and running the analysis (includes compiling of relevant __Stan__ code). 

	 usage: run_luxPom.py -t REGION_TOTAL_FILE -m REGION_METHYLATED_FILE -d DESIGN_MATRIX -o OUTFOLDER -a BS_EFF -b BS_BEFF -c SEQ_ERR -l $STAN_HOME
	 
	 Estimates methylation levels and infers differential methylation
	 
	 optional arguments:
	 -h, --help								show this help message and exit
	 -t REGION_TOTAL_FILE, --total_file REGION_TOTAL_FILE			file containing total read counts in regions; if not supplied, defaults to hmm_output/total_reads_all.txt
	 -m REGION_METHYLATED_FILE, --methylated_file REGION_METHYLATED_FILE	file containing methylated read counts in regions; if not supplied, defaults to hmm_output/total_reads_all.txt
	 -d DESIGN_MATRIX, --design_matrix DESIGN_MATRIX			file containing design matrix
	 -o OUTFOLDER, --outfolder OUTFOLDER					directory containing data analysis output; if not supplied, defaults to results
	 -a BS_EFF, --bs_Eff BS_EFF						comma-delimited bisulfite conversion efficiencies of replicates; if not supplied, defaults to 1.0
	 -b BS_BEFF, --bs_BEff BS_BEFF						comma-delimited incorrect bisulfite conversion efficiencies of replicates; if not supplied, defaults to 0
	 -c SEQ_ERR, --seq_Err SEQ_ERR						comma-delimited sequencing error rates of replicates; if not supplied, defaults to 0

	 -l $STAN_HOME, --cmdstan_loc $STAN_HOME 				CmdStan directory with full pathname
	 -v, --version								show program's version number and exit

For instance, run_luxPom.py can be called as

	python run_luxPom.py -t hmm_output/total_reads_all.txt -m hmm_output/methylated_reads_all.txt -d design_matrix.txt-l $STAN_HOME

*Input*

* **total\_reads\_all.txt** and **methylated\_reads\_all.txt** contain the combined total and methylated read counts, respectively, with one line per region. These files are returned by **run\_hmm.py** described above. If supplied, BS_EFF, BS_BEFF and SEQ_ERR must be in the same order as these files.
* **design\_matrixt.txt** contains the design matrix. The covariate to be tested must be in the third column (after rownames and intercept). Header contains the names of the covariates and row names contain the names of the replicates. The names of the replicates must be in the same order as the columns of **total\_reads\_all.txt** and **methylated\_reads\_all.txt** and, if supplied, BS_EFF, BS_BEFF and SEQ_ERR.

*Output*
* **regions.bed** is a **bed** file with 4 columns: 1) chromosome, 2) start position, 3) end position and 4) Bayes factor. Bayes factors are computed from the approximate posteriors of the coefficient of interest using the Savage-dickey method.


**References**

[1] J. Schreiber, “Pomegranate: fast and flexible probabilistic modeling in python.,” Journal of Machine Learning Research, 18.164:1-6, 2018.

[2] T. Äijö, X. Yue, A. Rao and H. Lähdesmäki, “LuxGLM: a probabilistic covariate model for quantification of DNA methylation modifications with complex experimental designs.,” Bioinformatics, 32.17:i511-i519, Sep 2016. 

[3] V. Halla-Aho & H. Lähdesmäki, “LuxUS: DNA methylation analysis using generalized linear mixed model with spatial correlation.,” Bioinformatics, 36.17:4535-4543, 2020.

[4] F. Krueger and S. R. Andrews, “Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications.,” Bioinformatics, 27.11:1571-1572, Jun 2011. 

[5] Stan Development Team, “Stan Modeling Language Users Guide and Reference Manual, Version 2.14.0.,” http://mc-stan.org, Jan 2016. 

[6] Stan Development Team, “PyStan: the Python interface to Stan, Version 2.14.0.0.,” http://mc-stan.org, Jan 2016.

