import os
import math
import glob
from pandas import DataFrame
import matplotlib.pyplot as plt
plt.switch_backend('agg')

colors = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"]

def merge_profiles(outdir, var, value=False):
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
		exit(1)

	if value:
		return df, prefix, out_name


def merge_barplot(outdir, var, tissue):
	
	df, prefix, out_name = merge_profiles(outdir, var, value=True)
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
	elif tissue == 'sperm':
		df = df.reindex(columns = ['miRNA','tsRNA','rsRNA','snoRNA','lincRNA','mRNA','others'])

	df_pcts = df.div(df.sum(1).astype(float), axis=0) * 100
	percent_out = prefix + out_name.split('.')[0] + '_percent.pdf'
	percent_txt = prefix + out_name.split('.')[0] + '_percent.txt'
	df_pcts.to_csv(percent_txt, header=True, sep='\t', float_format='%.2f')
	df_pcts.plot(kind='bar', fontsize=fontsize, width=0.8,stacked=True, color=colors, figsize = (width, height))
	plt.legend(loc = 0,borderaxespad = 0.5,fontsize = 'large',bbox_to_anchor = (1,1))
	plt.ylabel("Percent", fontsize=20)
	plt.savefig(percent_out, bbox_inches='tight')
	plt.close()


