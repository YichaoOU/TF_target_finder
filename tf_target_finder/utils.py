import argparse
import uuid
import pandas as pd
import os
from collections import OrderedDict
import sys
import numpy as np
import datetime
import getpass
import re, string
import warnings
import scipy
warnings.filterwarnings("ignore")


class query:
	def __init__(self,query_bed):
		df = pd.read_csv(query_bed,sep="\t",header=None)
		self.query_bed_Ncols = df.shape[1]
		df["query_chr"] = df[0]
		df["query_start"] = df[1].astype(int)
		df["query_end"] = df[2].astype(int)
		df['query_name'] = df[0]+":"+df["query_start"].astype(str)+"-"+df["query_end"].astype(str)
		self.query_bed = query_bed
		self.bed_cols = ['query_chr','query_start','query_end']
		self.df = df[self.bed_cols+['query_name']]
		self.log=OrderedDict()
		self.log['Number of query regions'] = df.shape[0]

	def extend_query(self,d):
		self.df["query_extend_start"] = self.df['query_start']-d
		self.df["query_extend_end"] = self.df['query_end']+d
	def find_nearest_TSS(self,TSS_bed):
		df = bedtools_closest(self.query_bed,TSS_bed)
		df.index = df[0]+":"+df[1].astype(str)+"-"+df[2].astype(str)
		self.df['nearest_TSS_gene'] = df.loc[self.df.query_name.tolist()][df.columns[self.query_bed_Ncols+4-1]].tolist()
		self.df['nearest_TSS_distance'] = df.loc[self.df.query_name][df.columns[-1]].tolist()
		
	def is_promoter(self,TSS_flank):
		# self.df['is_promoter'] = self.df['nearest_TSS_distance'] <= TSS_flank
		self.df['hard_assignment'] = self.df.apply(lambda x:get_TSS_gene(x,TSS_flank),axis=1)

	def find_EPI_target(self,EPI_bed,label,extend_EPI):
		df = pd.read_csv(EPI_bed,sep="\t",header=None)
		df[1] = df[1].astype(int)-extend_EPI
		df[2] = df[2].astype(int)+extend_EPI
		score_col_flag = df.shape[1]>=5
		df = bedtools_overlap(self.df[['query_chr','query_extend_start','query_extend_end','query_name']],df[df.columns[:5]])
		df.index = df[3].tolist()
		# print (df.head())
		use_col=-2
		if score_col_flag:
			use_col = -3
			
		target_genes = pd.DataFrame(df.groupby(3)[df.columns[use_col]].apply(lambda x: ','.join(x)))
		target_genes.columns = ['target_gene']
		target_genes = target_genes['target_gene'].to_dict()
		self.df['%s_gene'%(label)] = self.df['query_name'].map(target_genes)
		self.log['Number of assigned query regions by %s'%(label)] = len(target_genes)
		
		if score_col_flag:
			target_score = pd.DataFrame(df.groupby(3)[df.columns[-2]].apply(lambda x: ",".join([str(i) for i in x])))
			target_score.columns = ['target_score']
			target_score = target_score['target_score'].to_dict()		
			self.df['%s_score'%(label)] = self.df['query_name'].map(target_score)
		self.df= self.df.fillna(".")
	def get_final_target_genes(self,use_cols,label):
		self.df['all_target_genes_%s'%(label)] = self.df.apply(lambda x:combine_rows(x,use_cols),axis=1)
		out_file = "%s.query.targets_all.bed"%(label)
		self.df[self.bed_cols+['all_target_genes_%s'%(label)]].to_csv(out_file,index=False,header=False,sep="\t")
		
	def deg_evidence_filter(self,label,RNA_seq,FDR_col_name,LFC_col_name,FDR_cutoff, LFC_cutoff):
		"""add new column to subset the EPI target		
		deg filter		
		"""
		self.log['DEG FDR cutoff'] = FDR_cutoff
		self.log['DEG LFC cutoff'] = LFC_cutoff
		df = pd.read_csv(RNA_seq,sep="\t",index_col=0)
		if FDR_col_name != "FDR":
			df['FDR'] = df[FDR_col_name]
		if LFC_col_name != "LFC":
			df['LFC'] = df[LFC_col_name]
		self.deg_tsv = df
		self.FDR_dict = df[FDR_col_name].to_dict()
		self.LFC_dict = df[LFC_col_name].to_dict()
		self.df['targets_deg_filter_%s'%(label)] = self.df['all_target_genes_%s'%(label)].apply(lambda x:get_exp(self.FDR_dict,self.LFC_dict,FDR_cutoff, LFC_cutoff,x))
		self.df['Num_targets_deg_filter_%s'%(label)] = [len(x) for x in self.df['targets_deg_filter_%s'%(label)]]
		self.log['Number of assigned query regions by %s'%(label)] = self.df[self.df['Num_targets_deg_filter_%s'%(label)]>0].shape[0]
		out_file = "%s.query.DEG_targets_filter.bed"%(label)
		self.df[self.df['Num_targets_deg_filter_%s'%(label)]>0][self.bed_cols+['targets_deg_filter_%s'%(label)]].to_csv(out_file,index=False,header=False,sep="\t")
		
	def save_filtered_deg_table(self,label):
		myList = self.df['targets_deg_filter_%s'%(label)].tolist()
		flat_list = [item for sublist in myList for item in sublist]
		out_file = '%s.deg_table.tsv'%(label)
		self.deg_tsv.loc[flat_list].to_csv(out_file,sep="\t")	
		self.log['Number of assigned targets by %s'%(label)] = len(flat_list)
		
	def print_log(self):
		for k in self.log:
			print ("%s: %s"%(k,self.log[k]))
		
def get_TSS_gene(r,TSS_flank):
	if r.nearest_TSS_distance <= TSS_flank:
		return r.nearest_TSS_gene
	else:
		return "."
		
def combine_rows(r,use_cols):
	out = []
	for c in use_cols:
		# print (r.name,c,r[c])
		out+=r[c].split(",")
	out = list(set(out))
	try:
		out.remove(".")
	except:
		pass
	return out
			

def get_exp(FDR_dict,LFC_dict,FDR_cutoff, LFC_cutoff,x):
	# print (x)
	myList = x
	out = []
	for g in myList:
		try:
			FDR = FDR_dict[g]
			LFC = LFC_dict[g]
			if abs(LFC)>=LFC_cutoff and FDR<=FDR_cutoff:
				out.append(g)			
		except:
			continue
	return out
	

	
def bedtools_overlap(df1,df2):
	"""perform bedtools intersect for df1 and df2"""
	f1 = str(uuid.uuid4()).split("-")[-1]
	f2 = str(uuid.uuid4()).split("-")[-1]
	to_bed(df1,f1)
	to_bed(df2,f2)
	output = str(uuid.uuid4()).split("-")[-1]
	command = "bedtools intersect -a %s -b %s -wo > %s"%(f1,f2,output)
	os.system(command)
	
	df = pd.read_csv(output,sep="\t",header=None)
	# print (df.shape)
	delete_file(f1)
	delete_file(f2)
	delete_file(output)
	# print (df.head())
	return df	
def sort_bed(input,output):
	os.system("sort -k1,1 -k2,2n %s>%s"%(input,output))
		
def bedtools_closest(f1,f2):
	"""perform bedtools closest for input file 1 and 2"""
	outFile = str(uuid.uuid4()).split("-")[-1]
	sort1 = str(uuid.uuid4()).split("-")[-1]
	sort2 = str(uuid.uuid4()).split("-")[-1]
	sort_bed(f1,sort1)
	sort_bed(f2,sort2)
	command = "bedtools closest -a %s -b %s -d > %s"%(sort1,sort2,outFile)
	os.system(command)
	
	df = pd.read_csv(outFile,sep="\t",header=None)
	# print (df.head())
	delete_file(outFile)
	delete_file(sort1)
	delete_file(sort2)
	df['name'] = df[0]+":"+df[1].astype(str)+"-"+df[2].astype(str)
	df = df.drop_duplicates('name')	
	df = df.drop(['name'],axis=1)
	
	return df
	
	
def to_bed(df,out):
	df.to_csv(out,sep="\t",header=False,index=False)

def delete_file(f):
	os.system("rm %s"%(f))

def read_bed(f,d):
	df = pd.read_csv(f,sep="\t",header=None)
	df[1] = df[1].astype(int)
	df[2] = df[2].astype(int)
	df[1] = df[1]-d
	df[1][df[1] < 0] = 0
	df[2] = df[2]+d
	return df


def wccount(filename):
	df = pd.read_csv(filename,sep="\t",header=None)
	df['name'] = df[0]+":"+df[1].astype(str)+"-"+df[2].astype(str)
	df = df.drop_duplicates("name")
	# out = subprocess.Popen(['wc', '-l', filename],stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0]
	return df.shape[0]

