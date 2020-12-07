#!/usr/bin/env python3

### check python version
from sys import exit, version_info
if not version_info.major == 3 and version_info.minor >= 6:
	#print(version_info)
	print("\nThis script requires Python 3.6 or higher!")
	print("\nYou are using Python {}.{}\n".format(version_info.major, version_info.minor))
	exit(1)

from sys import argv, path
import os
import subprocess
import argparse
import datetime
import time
import glob
import re

from collections import defaultdict
from pandas import Series, DataFrame
from Bio import SeqIO

import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')

### import my library
from general import *
from common import *
from testis import *
from EV import *
from sperm import *
from merge import *
colors = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"]

###############
def srat(finput, outdir, library, tissue, threads, no_merge):
	
	### get the data with fastq format
	files=[]
	fAll = glob.glob("%s/*" % finput)
	for f in fAll:
		if f.endswith('fastq') or f.endswith('fastq.gz') or f.endswith('fq.gz') or f.endswith('fq'):
			files.append(f)
	
	if len(files) == 0 :
		print("\nNo fastq format files, please check your input files!\n")
		exit(1)
	
	### processing
	files=sorted(files)
	for i in files:
		dir_name = i.split("/")[-1].split("_")[0]
		directory = os.path.join(outdir, dir_name)

		### make directory
		if not os.path.exists(directory):
			try:
				os.makedirs(directory)
			except Exception as e:
				pass

		inputfile = i
		prefix = os.path.join(directory, dir_name)
		
		### cutadapt
		cut_adapt, collapser, trimmed = cutadapter(prefix, inputfile, threads)

		###
		devnull = open(os.devnull, 'w')
		### Common analysis
		if tissue == 'common':
			# 1.
			bowtie_spikein_out, bowtie_out_combined, map_genome, unmap_genome = common(prefix, library, threads,collapser, devnull)
			# 2.
			RNA_length, total_RNA, dic_miR, dic_miR_5p, dic_type, sum_spikein = commonProcess(bowtie_out_combined, prefix, cut_adapt, collapser, map_genome, bowtie_spikein_out)
			# 3.
			commonPlot(RNA_length,total_RNA, prefix)
			# 4.
			genomePlot(map_genome, unmap_genome, prefix)
			# 5.
			miRNAanalysis(prefix, dic_miR, dic_miR_5p, dic_type, sum_spikein, dir_name)
		
		### Testis
		elif tissue == 'testis':
			# 1.
			bowtie_spikein_out, bowtie_out_combined, map_genome, unmap_genome = testis(prefix, library, threads, collapser, devnull)
			# 2.
			RNA_length, total_RNA, dic_miR, dic_miR_5p, dic_type, sum_spikein = commonProcess(bowtie_out_combined, prefix, cut_adapt, collapser, map_genome, bowtie_spikein_out)
			# 3.
			testis_piRNA(bowtie_out_combined, prefix)
			# 4.
			testisPlot(RNA_length,total_RNA, prefix)
			# 5.
			genomePlot(map_genome, unmap_genome, prefix)
			# 6.
			miRNAanalysis(prefix, dic_miR, dic_miR_5p, dic_type, sum_spikein, dir_name)
		
		### EV
		elif tissue == 'EV':
			# 1.
			bowtie_spikein_out, bowtie_out_combined, map_genome, unmap_genome = EV(prefix, library, threads, collapser, devnull)
			# 2.
			RNA_length, total_RNA, dic_miR, dic_miR_5p, dic_type, sum_spikein = commonProcess(bowtie_out_combined, prefix, cut_adapt, collapser, map_genome, bowtie_spikein_out)
			# 3.
			dic_YRNA_type = EV_YRNA(bowtie_out_combined, prefix)
			# .
			EV_Plot(RNA_length,total_RNA, prefix, dic_YRNA_type)
			# 4.
			genomePlot(map_genome, unmap_genome, prefix)
			# 5.
			miRNAanalysis(prefix, dic_miR, dic_miR_5p, dic_type, sum_spikein, dir_name)
		
		elif tissue == 'sperm':
			# 1.
			bowtie_spikein_out, bowtie_out_combined, map_genome, unmap_genome = sperm(prefix, library, threads, collapser, devnull)
			# 2.
			sperm_RNA(bowtie_out_combined, prefix)
			# 3.
			RNA_length, total_RNA, dic_miR, dic_miR_5p, dic_type, sum_spikein = commonProcess(bowtie_out_combined, prefix, cut_adapt, collapser, map_genome, bowtie_spikein_out)
			# 4.
			spermPlot(RNA_length,total_RNA, prefix)
			# 5.
			genomePlot(map_genome, unmap_genome, prefix)
			# 6.
			miRNAanalysis(prefix, dic_miR, dic_miR_5p, dic_type, sum_spikein, dir_name)


		subprocess.call("rm %s %s" % (trimmed, cut_adapt), shell=True)
		devnull.close()
		print("\n%s ..... Finished %s" % (current_date(), dir_name))
	
	# ======================= #
	# merge expression files
	if no_merge:
		merge_profiles(outdir, "*summary_count.txt")
		merge_profiles(outdir, "*miRNA_counts_5p.txt")
		merge_profiles(outdir, "*miRNA_counts.txt")
		merge_barplot(outdir, "*pie_RNA.txt", tissue=args.tissue)
		
		if tissue == 'sperm':
			merge_profiles(outdir, "*tsRNA_counts.txt")
			merge_profiles(outdir, "*rsRNA_counts.txt")
			merge_profiles(outdir, "*piRNA_seq.txt")

		elif tissue == 'testis':
			merge_profiles(outdir, "*piRNA_seq.txt")

		elif tissue == "EV":
			merge_profiles(outdir, "*YRNA_counts.txt")
			
	# ======================= #
	# for spikein 


# =============================================== #

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'A small RNA analysis and visualization tool')

	parser.add_argument('-i', '--input', required=True, help='the input directory of raw data')
	parser.add_argument('-o', '--outdir', default='processed_data', help='output directory')
	parser.add_argument('-l', '--library', required=True, type=str, help='the reference sequence for mapping and annotation')
	parser.add_argument('-t', '--tissue', default='common', choices=['common', 'testis', 'sperm', 'EV'], type=str, help='the sample tissue type')
	parser.add_argument('-p', '--threads', default=1, type=int, help='number of alignment threads to launch (default: 1)')
	parser.add_argument('--no_merge', action='store_false', help="not merge the expression files")
	parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1')

	args = parser.parse_args()

	### Start
	print("\n%s ..... Start small RNA processing" % (current_date()))
	start_time = time.time()

	srat(args.input, args.outdir, args.library, args.tissue, args.threads, args.no_merge)
	
	end_time = time.time()
	run_time = round((end_time - start_time)/60, 5)
	### End
	print("\n%s ..... Finished all. Used time: %s m\n" % (current_date(), run_time))
