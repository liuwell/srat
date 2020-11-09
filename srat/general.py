
from sys import argv, path, exit
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


#########################
def cutadapter(prefix, inputfile):

	### cut adapter
	cut_adapt = prefix + ".cutadapt.fastq"
	cut_log = prefix + ".cutadapt.log"
	subprocess.call("cutadapt -a TGGAATTCTCGGGTGCCAAGG -O 10 -o %s %s > %s" % (cut_adapt, inputfile, cut_log), shell=True)

	### reads length distribution
	plotname = cut_adapt + ".png"
	writename = cut_adapt + ".txt"
	fastq_length_plot(cut_adapt, plotname, writename)

	### trimming short reads
	trimmed = prefix + ".tmp.fastq"
	collapser = prefix + ".collapser.fa"

	subprocess.call("cutadapt -q 15 -m 17 -f fastq -M 60 %s -o %s >> %s" % (cut_adapt, trimmed, cut_log), shell=True)

	subprocess.call("cat %s|fastq_to_fasta |fastx_collapser -o %s" % (trimmed, collapser), shell=True)
	return cut_adapt, collapser, trimmed


#########################
def genomePlot(map_genome, unmap_genome, prefix):	
	### length plot, only mapped genome
	d_series = Series(get_fasta_length(map_genome))
	plot_map = prefix + ".length_mapGenome.png"
	df = d_series.sort_index()
	df.plot(kind = "bar", color = "#990000", fontsize = 15, width = 1)
	plt.xticks(range(0,35,5),("17", "22", "27", "32", "37", "42", "47"), fontsize=15, rotation=0)
	plt.xlabel("Length",fontsize=15)
	plt.ylabel("Counts",fontsize=15)
	plt.savefig(plot_map, bbox_inches='tight')
	plt.close()

	### length plot, unmapped genome
	d_series = Series(get_fasta_length(unmap_genome))
	plot_unmap = prefix + ".length_unmapGenome.png"
	df = d_series.sort_index()
	df.plot(kind = "bar", color = "#990000", fontsize = 15, width = 1)
	plt.xticks(range(0,35,5),("17", "22", "27", "32", "37", "42", "47"), fontsize=15, rotation=0)
	plt.xlabel("Length",fontsize=15)
	plt.ylabel("Counts",fontsize=15)
	plt.savefig(plot_unmap, bbox_inches='tight')
	plt.close()


##########################
def miRNAanalysis(prefix, dic_miR, dic_miR_5p, dic_type, sum_spikein, dir_name):
	### miRNA counts
	mir_out = open(prefix + ".miRNA_counts.txt", 'w')
	dic_miR = sorted(dic_miR.items(), key=lambda d:d[1], reverse=True)
	for k in dic_miR:
		mir_out.write(k[0] + "\t"+ str(k[1]) + "\n")
	mir_out.close()
	
	### normalized by total mapped reads
	total_map = dic_type["total_map"]
	mir_file2 = prefix + ".miRNA_mapped.txt"
	mir_out2 = open(mir_file2,'w')
	for k in dic_miR:
		mir_out2.write(k[0] + "\t"+ str(int(k[1]/total_map * 1000000)) + "\n")
	mir_out2.close()

	### normalized by total miRNAs
	total_miR = dic_type["miRNA"]
	mir_file3 = prefix + ".miRNA_miR.txt"
	mir_out3 = open(mir_file3,'w')
	for k in dic_miR:
		mir_out3.write(k[0] + "\t"+ str(int(k[1]/total_miR * 1000000)) + "\n")
	mir_out3.close()

	### miRNA counts, no 5p isoform
	mir_out4 = open(prefix+".miRNA_counts_5p.txt", 'w')
	dic_miR_5p = sorted(dic_miR_5p.items(), key=lambda d:d[1], reverse=True)
	for k in dic_miR_5p:
		mir_out4.write(k[0] + "\t"+ str(k[1]) + "\n")
	mir_out4.close()

	
	### nomalized by spikein
	if sum_spikein > 0:
		mir_spike = prefix + ".miRNA_spikein.txt"
		mir_spike_out = open(mir_spike, 'w')
		for k in dic_miR:
			mir_spike_out.write(k[0] + '\t' + str(int(k[1]/sum_spikein * 100000)) + '\n')
		mir_spike_out.close()
	else:
		print("\n%s ..... No spikein %s" % (current_date(), dir_name))
	


################
def current_date():
	return datetime.datetime.now().strftime('%b-%d-%Y %H:%M:%S')

################
### get total count of fasta format file
def get_count(data):
	with open(data) as handle:
		reads = 0
		lock =re.search("fastq",data)
		for line in handle:
			if lock:
				reads += 0.25
			else:
				if re.match(">", line):
					counts = line.split("-")[-1]
					reads += int(counts)
		return int(reads)

################
### get the length distribution of fasta format file
def get_fasta_length(data):
	with open(data) as handle:
		d = defaultdict(int)
		for line in SeqIO.parse(handle, 'fasta'):
			d[len(line.seq)] += int(line.id.split("-")[1])
		return d

################
### fastq length plot
def fastq_length_plot(fastq, plotname, writename):
	dict_length = defaultdict(int)
	with open(fastq) as handle:
		for record in SeqIO.parse(handle, "fastq"):
			l = len(record.seq)
			dict_length[l] += 1

		for i in range(0,150):
			if i in dict_length.keys():
				pass
			else:
				dict_length[i] = 0

	df = Series(dict_length)
	df = df.sort_index()
	#print df
	df.to_csv(writename, header=False, sep="\t", float_format="%.0f")
	df.plot(kind = "bar", color = "#990000", fontsize = 15, width = 1)
	plt.xlim(-1, 61)
	plt.xticks(range(0, 61, 10), ("0", "10", "20", "30", "40", "50", "60"), rotation=0)
	plt.savefig(plotname, bbox_inches='tight')
	plt.close()

