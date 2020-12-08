#!/usr/bin/env python3

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

from general import *

colors = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"]

def EV(prefix, library, threads, collapser, devnull):

	# for spikein	
	#bowtie_spikein_out = prefix + ".bowtie_spikein.out"
	#unmapped_spikein = prefix + ".unmap_spikein.fa"
	#subprocess.call("bowtie -v 0 %s/../spikein/index -p %d -f %s --norc --un %s > %s" % (library, threads, collapser, unmapped_spikein, bowtie_spikein_out), shell=True, stderr=devnull)
	
	# for miRNA
	unmapped_mir = prefix + ".unmap_mir.fa"
	bowtie_mir_out = prefix + ".bowtie_miRNA.out"
	subprocess.call("bowtie -v 0 %s/index_miRNA --norc --suppress 6,7 -p %d -f %s --un %s > %s" % (library, threads, collapser, unmapped_mir, bowtie_mir_out), shell=True, stderr=devnull)
	# for YRNA
	unmapped_YRNA = prefix + ".unmap_YRNA.fa"
	bowtie_YRNA_out = prefix + ".bowtie_YRNA.out"
	subprocess.call("bowtie -v 0 %s/index_YRNA --norc --suppress 6,7 -p %d -f %s --un %s > %s" % (library, threads, unmapped_mir, unmapped_YRNA, bowtie_YRNA_out), shell=True, stderr=devnull)
	# for other RNA
	unmapped_RNA = prefix + ".unmap_RNA.fa"
	bowtie_RNA_out = prefix + ".bowtie_RNA.out"
	subprocess.call("bowtie -v 0 %s/index_RNA --norc --suppress 6,7 -p %d -f %s --un %s > %s" % (library, threads, unmapped_YRNA, unmapped_RNA, bowtie_RNA_out), shell=True, stderr=devnull)
	# for genome
	map_genome = prefix + ".map_genome.fa"
	unmap_genome = prefix + ".unmap_genome.fa"
	subprocess.call("bowtie -v 0 %s/genome -p %d -f %s --al %s --un %s " % (library, threads, unmapped_RNA, map_genome, unmap_genome), shell=True, stderr=devnull, stdout=devnull)
	# combined
	bowtie_out_combined = prefix + ".bowtie_combined.out"
	os.system("cat %s %s %s> %s" % (bowtie_mir_out, bowtie_YRNA_out, bowtie_RNA_out, bowtie_out_combined))
	os.system("rm %s %s %s" % (unmapped_mir, unmapped_YRNA,unmapped_RNA))
	os.system("rm %s %s %s" % (bowtie_mir_out, bowtie_YRNA_out, bowtie_RNA_out))

	return bowtie_out_combined, map_genome, unmap_genome

# ========================== #
def EV_YRNA(bowtie_out_combined, prefix):
	
	dic_YRNA = defaultdict(int)
	dic_YRNA_type = defaultdict(int)
	with open(bowtie_out_combined) as handle:
		for line in handle:
			seg =line.split()
			count = int(seg[0].split("-")[1])
			
			if re.search("YRNA", seg[2]):
				yR = seg[4] + "|" + seg[2]
				dic_YRNA[yR] += count
				yR_type = seg[2].split("|")[1]
				dic_YRNA_type[yR_type] += count

	YRNA_out = prefix + ".YRNA_counts.txt"
	YRNA_series = Series(dic_YRNA)
	YRNA_series.to_csv(YRNA_out, header=False, sep='\t')
	return dic_YRNA_type

# ========================== #
def EV_Plot(dic_length_RNA, total_RNA, prefix, dic_YRNA_type):

	### length distribution of different RNA types
	#df = DataFrame(dic_length_RNA).T
	df = dic_length_RNA.fillna(value=0)
	df = df.sort_index(axis=0, ascending=True)
	df_RNA = df.loc[:,['miRNA', 'YRNA', 'tsRNA', 'rsRNA', 'snoRNA', 'lncRNA', 'mRNA']]

	df_RNA_other = pd.DataFrame(df.sum(axis=1) - df_RNA.sum(axis=1),columns=["others"])
	df_RNA = df_RNA.join(df_RNA_other)

	plot_RNA = prefix + ".length_RNA_counts.pdf"
	df_RNA.plot(kind='bar', stacked=True, fontsize = 15, color=colors, width=1, linewidth=0.01)
	plt.xticks(range(0,35,5),("17", "22", "27", "32", "37", "42", "47"), fontsize=15, rotation=0)
	plt.xlabel("Length",fontsize=15)
	plt.ylabel("Counts",fontsize=15)
	plt.savefig(plot_RNA, bbox_inches='tight')
	plt.close()
	df_RNA_csv = prefix + ".length_RNA_counts.txt"
	df_RNA.to_csv(df_RNA_csv, header=False, sep='\t')
	
	### length plot percent
	df_RNA_percent = df_RNA/total_RNA * 100
	plot_RNA_percent = prefix + ".length_RNA_percent.pdf"
	df_RNA_percent.plot(kind='bar', stacked=True, figsize=(5.5,4), fontsize = 15, color=colors, width=1, linewidth=0.01)
	plt.xticks(range(0,35,5),("17", "22", "27", "32", "37", "42", "47"), fontsize=15, rotation=0)
	plt.xlabel("Length",fontsize=15)
	plt.ylabel("Percent",fontsize=15)
	plt.savefig(plot_RNA_percent, bbox_inches='tight')
	plt.close()

	### pie plot
	pie_RNA = prefix + ".pie_RNA.pdf"
	df_RNA_sum = df_RNA.sum(axis=0)
	df_RNA_sum.name = ''
	df_RNA_sum.plot(kind='pie', figsize=(6,6), colors=colors, autopct='%.1f',fontsize=15 )
	pie_csv = prefix + ".pie_RNA.txt"
	df_RNA_sum = df_RNA_sum.fillna(value=0)
	df_RNA_sum.to_csv(pie_csv, header=False, sep='\t')
	plt.savefig(pie_RNA, bbox_inches='tight')
	plt.close()


	### pie plot of YRNA type
	df_YRNA = Series(dic_YRNA_type)
	df_YRNA = df_YRNA.fillna(value=0)
	df_YRNA = df_YRNA.reindex(['Y5', 'Y4', 'Y3', 'Y1'], fill_value=0)
	pie_YRNA = prefix + ".pie_YRNA.pdf"
	df_YRNA.name = ''
	df_YRNA.plot(kind='pie', figsize=(4,4), colors=colors, autopct='%.1f',fontsize=15 )
	plt.savefig(pie_YRNA, bbox_inches='tight')
	plt.close()
	
	pie_YRNA_csv = prefix + ".pie_YRNA.txt"
	df_YRNA.to_csv(pie_YRNA_csv, header=False, sep='\t')


