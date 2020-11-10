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

colors = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"]

###############
def srat(finput, outdir, library, tissue, threads):
	
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

		### path
		#pre_index = path[0] + '/../smallRNA_index' ##### set the position of RNA index file
		#pre_index = library
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
			RNA_length, total_RNA, dic_miR, dic_miR_5p, dic_type, sum_spikein = testisProcess(bowtie_out_combined, prefix, cut_adapt, collapser, map_genome, bowtie_spikein_out)
			# 3.
			testisPlot(RNA_length,total_RNA, prefix)
			# 4.
			genomePlot(map_genome, unmap_genome, prefix)
			# 5.
			miRNAanalysis(prefix, dic_miR, dic_miR_5p, dic_type, sum_spikein, dir_name)
		
		### EV
		elif tissue == 'EV':
			# 1.
			bowtie_spikein_out, bowtie_out_combined, map_genome, unmap_genome = EV(prefix, library, threads, collapser, devnull)
			# 2.
			RNA_length, total_RNA, dic_miR, dic_miR_5p, dic_type, sum_spikein = EV_Process(bowtie_out_combined, prefix, cut_adapt, collapser, map_genome, bowtie_spikein_out)
			# 3.
			EV_Plot(RNA_length,total_RNA, prefix)
			# 4.
			genomePlot(map_genome, unmap_genome, prefix)
			# 5.
			miRNAanalysis(prefix, dic_miR, dic_miR_5p, dic_type, sum_spikein, dir_name)
		
		subprocess.call("rm %s %s" % (trimmed, cut_adapt), shell=True)
		devnull.close()
		print("\n%s ..... Finished %s" % (current_date(), dir_name))



############################
### merge RNA profiles
import math
def merge_profiles(outdir, var, barplot=False, tissue='common'):
	files = glob.glob("%s/*/%s" % (outdir,var))
	#print(files)
	if files  :

		files = sorted(files)
		out_name = files[0].split(".",1)[1]
		dict_merge = {}
		for f in files:
			with open(f) as handle:
				for line in handle:
					k_map = line.split()[0]
					k_RNA = f.split("/")[-2]
					count =float(line.split()[1])
					if k_map in dict_merge:
						dict_merge[k_map][k_RNA] = count
	
					else:
						tmp_dic = {}
						tmp_dic[k_RNA] = count
						dict_merge[k_map] = tmp_dic
		df = DataFrame(dict_merge).T
		df = df.fillna(value=0) ### fill NA to 0
		df_sum = DataFrame(df.sum(axis=1), columns=['sum'])
		df = df.join(df_sum)
		df = df.sort_values(by="sum", ascending=False) ### sort by sum
		df.drop(['sum'], axis=1, inplace=True)
	
		prefix = os.path.join(outdir, "merge_")
		merge_out = prefix + out_name
		df.to_csv(merge_out, sep="\t", header=True, float_format="%.0f")

	else:
		print("\n### Merge files %s failed, it is not exsit in %s/*/* \n" %(var, outdir))
		#exit(1)
	
	if barplot:
		df=df.T
		### ncol: df.shape[1], nrow: df.shape[0]
		width = int(df.shape[0])
		height = 6
		fontsize =20
		if width >= 8 :
			width = math.log(width, 2) * 2 ### adjust the width of barplot

		if tissue == 'common':
			df = df.reindex(columns = ['miRNA','tsRNA','rsRNA','snoRNA','misc_RNA','lincRNA','mRNA','others'])
		elif tissue == 'testis':
			df = df.reindex(columns = ['miRNA','piRNA','tsRNA','rsRNA','snoRNA','lincRNA','mRNA','others'])
		elif tissue == 'EV':
			df = df.reindex(columns = ['miRNA','YRNA','tsRNA','rsRNA','snoRNA','lincRNA','mRNA','others'])
		
		df_pcts = df.div(df.sum(1).astype(float), axis=0) * 100
		percent_out = prefix + out_name.split('.')[0] + '_percent.pdf'
		percent_txt = prefix + out_name.split('.')[0] + '_percent.txt'
		df_pcts.to_csv(percent_txt, header=True, sep='\t', float_format='%.2f')
		df_pcts.plot(kind='bar', fontsize=fontsize, width=0.8,stacked=True, color=colors, figsize = (width, height))
		plt.legend(loc = 0,borderaxespad = 0.5,fontsize = 'large',bbox_to_anchor = (1,1))
		plt.ylabel("Percent", fontsize=20)
		plt.savefig(percent_out, bbox_inches='tight')
		plt.close()




# =============================================== #

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'A small RNA analysis and visualization tool')

	parser.add_argument('-i', '--input', required=True, help='the input directory of raw data')
	parser.add_argument('-o', '--outdir', default='processed_data', help='output directory')
	parser.add_argument('-l', '--library', required=True, type=str, help='the reference sequence for mapping and annotation')
	#parser.add_argument('-s', '--species', required=True, choices=['human', 'mouse'], type=str, help='the small RNA reference')
	parser.add_argument('-t', '--tissue', default='common', choices=['common', 'testis', 'EV'], type=str, help='the sample tissue type')
	parser.add_argument('-p', '--threads', default=1, type=int, help='number of alignment threads to launch (default: 1)')

	#parser.add_argument('-piRNA', help = 'considering piRNA(mostly for germline cell)', action = 'store_true')
	#parser.add_argument('-YRNA', help = 'considering YRNA(mostly for human plasma EV)', action = 'store_true')

	args = parser.parse_args()

	### Start
	print("\n%s ..... Start small RNA processing" % (current_date()))
	start_time = time.time()

	srat(args.input, args.outdir, args.library, args.tissue, args.threads)
	
	#merge_profiles(args.outdir, "*miRNA_counts.txt")
	#merge_profiles(args.outdir, "*miRNA_counts_5p.txt")
	merge_profiles(args.outdir, "*miRNA_miR.txt")
	#merge_profiles(args.outdir, "*miRNA_mapped.txt")
	merge_profiles(args.outdir, "*summary_count.txt")
	merge_profiles(args.outdir, "*pie_RNA.txt", barplot=True, tissue=args.tissue)

	end_time = time.time()
	run_time = round((end_time - start_time)/60, 5)
	### End
	print("\n%s ..... Finished all. Used time: %s m\n" % (current_date(), run_time))
