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

def EV(prefix, pre_index, species, collapser, devnull):

	# for spikein	
	bowtie_spikein_out = prefix + ".bowtie_spikein.out"
	unmapped_spikein = prefix + ".unmap_spikein.fa"
	#index_spikein = '/home/liuwei/genome/smallRNA_index/spikein/index'
	subprocess.call("bowtie -v 0 %s/spikein/index -f %s --norc --un %s > %s" % (pre_index, collapser, unmapped_spikein, bowtie_spikein_out), shell=True, stderr=devnull)
	
	# for miRNA
	unmapped_mir = prefix + ".unmap_mir.fa"
	bowtie_mir_out = prefix + ".bowtie_miRNA.out"
	subprocess.call("bowtie -v 0 %s/%s/index_miRNA --norc --suppress 6,7 -f %s --un %s > %s" % (pre_index, species, unmapped_spikein, unmapped_mir, bowtie_mir_out), shell=True, stderr=devnull)
	# for YRNA
	unmapped_YRNA = prefix + ".unmap_YRNA.fa"
	bowtie_YRNA_out = prefix + ".bowtie_YRNA.out"
	subprocess.call("bowtie -v 0 %s/%s/index_YRNA --norc --suppress 6,7 -f %s --un %s > %s" % (pre_index, species, unmapped_mir, unmapped_YRNA, bowtie_YRNA_out), shell=True, stderr=devnull)
	# for other RNA
	unmapped_RNA = prefix + ".unmap_RNA.fa"
	bowtie_RNA_out = prefix + ".bowtie_RNA.out"
	subprocess.call("bowtie -v 0 %s/%s/index_RNA --norc --suppress 6,7 -f %s --un %s > %s" % (pre_index, species, unmapped_YRNA, unmapped_RNA, bowtie_RNA_out), shell=True, stderr=devnull)
	# for genome
	map_genome = prefix + ".map_genome.fa"
	unmap_genome = prefix + ".unmap_genome.fa"
	subprocess.call("bowtie -v 0 %s/%s/genome -f %s --al %s --un %s " % (pre_index, species, unmapped_RNA, map_genome, unmap_genome), shell=True, stderr=devnull, stdout=devnull)
	# combined
	bowtie_out_combined = prefix + ".bowtie_combined.out"
	os.system("cat %s %s %s > %s" % (bowtie_mir_out, bowtie_YRNA_out, bowtie_RNA_out, bowtie_out_combined))
	os.system("rm %s %s %s %s" % (unmapped_spikein, unmapped_mir, unmapped_YRNA, unmapped_RNA))

	return bowtie_spikein_out, bowtie_out_combined, map_genome, unmap_genome

#########################
def EV_Process(bowtie_out_combined, prefix, cut_adapt, collapser, map_genome, bowtie_spikein_out):

	dic_miR = defaultdict(int)         ### count miRNA 
	dic_miR_5p = defaultdict(int)      ### count miRNA, 5p in align
	dic_type = defaultdict(int)        ### count RNA type
	dic_all = defaultdict(int)         ### count RNA count
	dic_all_seq = defaultdict(int)     ### count RNA seq
	
	total_RNA = 0
	miRNA_kinds = set()

	dic_length_RNA = {}  ### length distribution of all different kinds mapped RNAs

	with open(bowtie_out_combined) as handle:
		for line in handle:
			seg =line.split()
			count = int(seg[0].split("-")[1])

			### miRNA
			if re.search("miRNA", seg[2]):
				mir = seg[2].split("|")[0]
				dic_miR[mir]+=count
				miRNA_kinds.add(mir)
				### no 5p isoform
				if seg[3] == "0":
					dic_miR_5p[mir]+=count
			
			title_split = seg[2].split("|")
			RNA_type = title_split[-1]
			dic_type[RNA_type] += count
			total_RNA += count

			k_map = len(seg[4])
			k_RNA = RNA_type
			if k_map in dic_length_RNA:
				if k_RNA in dic_length_RNA[k_map]:
					dic_length_RNA[k_map][k_RNA] += count

				else:
					dic_length_RNA[k_map][k_RNA] = count

			else:
				tmp_dic = {}
				tmp_dic[k_RNA] = count
				dic_length_RNA[k_map] = tmp_dic

			### piRNA
			if re.search("piRNA", seg[2]):
				pir_c = seg[2].split("|")[0]
				piRNA_cluster.add(pir_c)
				piRNA_seq.add(seg[4])
				dic_piR[seg[4]]+=count
				dic_piR_cluster[pir_c]+=count
	
	### summary RNA counts
	summary_out = open(prefix+'.summary_count.txt','w')
	dic_type["totals"] = get_count(cut_adapt)
	dic_type["useful"] = get_count(collapser)
	useful_ratio = round(dic_type["useful"]/dic_type["totals"]*100, 2)

	dic_type["total_map"] = get_count(map_genome)+total_RNA
	dic_type["total_knownRNA"] = total_RNA
	map_ratio = round(dic_type["total_map"]/dic_type["useful"]*100, 2)
	RNA_ratio = round(dic_type["total_knownRNA"]/dic_type["total_map"]*100, 2)
	dic_type["miRNA_kinds"] = len(miRNA_kinds)

	dic_type_sort = sorted(dic_type.items(), key=lambda d:d[1], reverse=True)
	for k in dic_type_sort:
		summary_out.write(k[0] + "\t"+ str(k[1]) + "\n")

	### for spikein
	with open(bowtie_spikein_out) as handle:
		sum_spikein = sum([int(i.split()[0].split("-")[1]) for i in handle])
	dic_type["spikein"] = sum_spikein
	spikein_ratio = round(dic_type["spikein"]/dic_type["total_map"]*100, 2)

	summary_out.write("useful_ratio(%)" + "\t" + str(useful_ratio) + "\n")
	summary_out.write("map_ratio(%)" + "\t" + str(map_ratio) + "\n")
	summary_out.write("knownRNA_ratio(%)" + "\t" + str(RNA_ratio) + "\n")
	summary_out.write("spikein" + "\t" + str(sum_spikein) + "\n")
	summary_out.write("spikein_ratio(map%)" + "\t" + str(spikein_ratio) + "\n")
	summary_out.close()
	
	### rename dic_length_RNA
	dic_length_RNA = DataFrame(dic_length_RNA).T
	dic_length_RNA.rename(columns={'tRNA':'tsRNA', 'rRNA':'rsRNA', 'protein_coding':'mRNA'}, inplace=True)

	return dic_length_RNA, total_RNA, dic_miR, dic_miR_5p, dic_type, sum_spikein


###################
def EV_Plot(dic_length_RNA, total_RNA, prefix):

	### length distribution of different RNA types
	#df = DataFrame(dic_length_RNA).T
	df = dic_length_RNA.fillna(value=0)
	df_RNA = df.loc[:,['miRNA', 'YRNA', 'tsRNA', 'rsRNA', 'snoRNA', 'lincRNA', 'mRNA']]

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
	#df_RNA_sum.name = dir_name
	df_RNA_sum.plot(kind='pie', figsize=(6,6), colors=colors, autopct='%.1f',fontsize=15 )
	pie_csv = prefix + ".pie_RNA.txt"
	df_RNA_sum = df_RNA_sum.fillna(value=0)
	df_RNA_sum.to_csv(pie_csv, header=False, sep='\t')
	plt.savefig(pie_RNA, bbox_inches='tight')
	plt.close()
