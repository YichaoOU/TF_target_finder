#!/usr/bin/env python


"""

Input
-----

1. gene expression table (3 columns)

2. OR known interaction table (2 columns)

3. Candidate gene list (1 column)

Output
-------

gene list




"""


import sys
import os
p_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
sys.path.append(os.path.abspath(p_dir+"../utils/"))
from utils import *



def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	mainParser.add_argument('-g',"--gene_list",help="1 column",required=True)	
	mainParser.add_argument('-e',"--evidence",help="1 or 2 or 3 columns, no header,1 = filter by given gene list, 2 = known interaction, 3 = gene expression",required=True)	
	mainParser.add_argument('-t',"--target",help="if input evidence is known interaction, you have to have a target gene, so that every interacting genes to this targets will be used as filter",default="None")	
	mainParser.add_argument("--cutoff",help="if input evidence is gene expression, please give a cutoff for logFC and FDR/p-value, e.g., 2,1e-3",default="2,1e-2")	
	mainParser.add_argument("--cols",help="columns of gene name, logFC and FDR/p-value, should give in order of name, logFC, FDR",default="None")	
	
	mainParser.add_argument('-o',"--output",  help="output overlapped gene list, [x].candidate.list",default="output")
	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args

def read_table(f):
	return pd.read_csv(f,sep="\t",header=None)

def evidence_to_list(args):
	df = read_table(args.evidence)
	if not args.cols == "None":
		df = pd.read_csv(args.evidence,sep="\t")
		df = df[args.cols.split(",")]
		print (df.head())
		df.columns = [0,1,2]
	if df.shape[1] == 1:
		return df[0].tolist()
	if df.shape[1] == 3:
		logFC,FDR = [float(x) for x in args.cutoff.split(",")]
		print ("using logFC cutoff of |%s|, and FDR/p-value cutoff of %s"%(logFC,FDR))
		df = df[df[1].abs() >= logFC]
		df = df[df[2] <= FDR]
		return df[0].tolist()
		
	if df.shape[1] == 2:
		df = df.dropna()
		df1 = df[df[0].str.contains(args.target,case=False)]
		df2 = df[df[1].str.contains(args.target,case=False)]
		df = pd.concat([df1,df2])
		return list(set(df[0]+df[1]))
	

def main():

	args = my_args()
	gene = read_table(args.gene_list)
	gene[0] = [x.upper() for x in gene[0]]
	evidence = evidence_to_list(args)
	evidence = [x.upper() for x in evidence]
	gene = gene[gene[0].isin(evidence)]
	print ("Number of genes remained: %s"%(gene.shape[0]))
	to_bed(gene,"%s.list"%(args.output))
	print ("Output file: %s.list"%(args.output))

if __name__ == "__main__":
	main()

































