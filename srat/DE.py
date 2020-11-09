#!/usr/bin/env python3
# python3.6
# ref link: https://www.jianshu.com/p/91c98585b79b
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import argparse

def DE(fi, wt, ko):

	prefix = fi.split('.')[0]

	data = pd.read_table(fi, header=0, index_col=0)
	data = np.log2(data+1)
	
	# Boxplot of the expression data
	color = {'boxes': 'DarkGreen', 'whiskers': 'DarkOrange', 'medians': 'DarkBlue', 'caps': 'Gray'}
	data.plot(kind='box', color=color, sym='r.', title=prefix)
	plt.xticks(rotation=40)
	plt.xlabel('Samples', fontsize=15)
	plt.ylabel('log2(CPM+1)', fontsize=15)

	out_box = prefix + "_boxplot.pdf"
	plt.savefig(out_box, bbox_inches='tight')
	plt.close()
	
	# Density plot of the expression data
	data.plot(kind='density', title=prefix)
	#data.plot.kde()
	out_density = prefix + "_density.pdf"
	plt.savefig(out_density)
	plt.close()
	
	
	#####
	wt = wt.split(':')
	wt1 = int(wt[0])
	wt2 = int(wt[1])
	ko = ko.split(':')
	ko1 = int(ko[0])
	ko2 = int(ko[1])

	# The mean expression of wt samples for each genes
	wt = data.iloc[:, wt1:wt2].mean(axis=1)
	# The mean expression of ko samples for each genes
	ko = data.iloc[:, ko1:ko2].mean(axis=1)
	
	# FoldChange, ko vs wt
	foldchange = ko - wt
	
	# P value
	pvalue = []
	gene_number = len(data.index)
	for i in range(0, gene_number):
		ttest = stats.ttest_ind(data.iloc[i,wt1:wt2], data.iloc[i, ko1:ko2])
		pvalue.append(ttest[1])
	
	
	### vocano plot
	pvalue_arr = np.asarray(pvalue)
	result = pd.DataFrame({'pvalue': pvalue_arr, 'FoldChange': foldchange})
	result['log10(pvalue)'] = -np.log10(result['pvalue'])
	
	result['sig'] = 'normal'
	result.loc[(result.FoldChange > 1)&(result.pvalue < 0.05), 'sig'] = 'up'
	result.loc[(result.FoldChange < -1)&(result.pvalue < 0.05), 'sig'] = 'down'
	
	sns.scatterplot(x="FoldChange", y="log10(pvalue)", hue='sig', hue_order=('down', 'normal', 'up'), palette=("#377EB8", "grey", "#E41A1C"), data=result)
	
	plt.xlabel('FoldChange')
	plt.ylabel('-log10(P value)')
	
	out_volcano = prefix + "_DE_volcano.pdf"
	plt.savefig(out_volcano)
	plt.close()

	out_result = prefix + "_FoldChange.txt"
	result.to_csv(out_result, header=True, sep='\t')

	### Heatmap
	fold_cutoff = 1
	pvalue_cutoff = 0.05
	
	filtered_ids=[]
	for i in range(0, gene_number):
		if(abs(foldchange[i]) >= fold_cutoff) and (pvalue[i] <= pvalue_cutoff):
			filtered_ids.append(i)
	
	filtered = data.iloc[filtered_ids, :]
	
	ncols = int(data.shape[1])
	#print(ncols)
	sns.clustermap(filtered, cmap='RdBu_r', standard_scale=0, figsize=(ncols, 10))
	
	out_heatmap = prefix + "_DE_heatmap.pdf"
	plt.savefig(out_heatmap)
	plt.close()

	print("\n# Finished the different expressed analysis")
	print("\n# Input file: %s" % fi)
	print("\n# Output txt: %s" % out_result)
	print("\n# Output figure: %s, %s, %s, %s\n" % (out_box, out_density, out_volcano, out_heatmap))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='For caculating differently expressed data, such as miRNA expression, gene expression')

	parser.add_argument('-i', '--input', required=True, help='the input data')
	parser.add_argument('-wt', required=True, type=str, help='the position of wildtype samples, such as 0:3')
	parser.add_argument('-ko', required=True, type=str, help='the position of knowkout samples, such as 3:6')

	args = parser.parse_args()

	DE(args.input, args.wt, args.ko)





