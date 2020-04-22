#!/usr/bin/env python

from tf_target_finder.utils import *

"""Assign targets given TSS annotation and promoter enhacer interaction

Input
-----

1. query set (bed 3, required)

This is a pandas dataframe object, all the subsequent analysis is just 

2. TSS annotation (bed 4, required)

3. Promoter enhancer interaction (bed 5, chr,start,end, gene, score)

4. Differentially expressed genes (gene name, FDR, logFC)

Output
-----

1. .query.targets_all.bed

2. .query.DEG_targets_filter.bed


"""

#!/usr/bin/env python


"""

This program uses bedtools to extract overlaped genes

Dependency
----------

Bedtools

Input
-----

TFBS bed file
Gene-associated features bed file

Output
-----

overlaped gene list

Parameters
---------

-d1: distance cutoff for TF bed

-d2: distance cutoff for Gene-associated features bed file

"""

def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	

	group = mainParser.add_mutually_exclusive_group(required=True)
	group.add_argument('-q',"--query_bed",  help="3 column bed file, additional columns are OK, but will be ignored")
	group.add_argument("-c","--config",  help="the config file is a tsv file containing all inputs and parameters")	

	mainParser.add_argument('-tss',"--tss_bed",help="4 column bed file, the 4th column should be gene name, should match to the gene name in DEG file (if supplied). Additional columns are OK, but will be ignored")		
	mainParser.add_argument("-epi","--epi_bed",help="5 column bed file, the 4th column should be gene name, should match to the gene name in DEG file and TSS annotation(if supplied). The 5th column should be score (optional). Additional columns are OK, but will be ignored")	
	
	mainParser.add_argument("-exp","--deg_tsv",help="any number of columns, first column should be gene name, first row should be column names. should contain FDR and LFC.")	
	
	mainParser.add_argument('--conservative',  help="only use captureC data for assignment", action='store_true')
	
	mainParser.add_argument("--LFC_col_name",help="LFC_col_name",default="logFC")	
	mainParser.add_argument("--FDR_col_name",help="FDR_col_name",default="adj.P.Val")	
	mainParser.add_argument("--LFC_cutoff",help="LFC cutoff",default=1,type=float)	
	mainParser.add_argument("--FDR_cutoff",help="FDR cutoff",default=0.05,type=float)	
	
	mainParser.add_argument("-d1",help="extend query bed for intersection",default=0,type=int)	
	mainParser.add_argument("-d2",help="extending tss for intersection",default=10000,type=int)	
	mainParser.add_argument("-d3",help="extending epi for intersection",default=5000,type=int)	
	

	mainParser.add_argument('-o',"--output",  help="output intermediate file",default="output")
	mainParser.add_argument("--label",  help="prefix for the file",default="TFBS")
	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args

	
def main():

	args = my_args()
	
	my_query = query(args.query_bed)
	my_query.extend_query(args.d1)
	my_query.find_nearest_TSS(args.tss_bed)
	my_query.is_promoter(args.d2)
	
	EPI_target_label = "captureC"
	my_query.find_EPI_target(args.epi_bed,EPI_target_label,args.d3)
	target_gene_set_label = args.label
	use_cols = ['nearest_TSS_gene','hard_assignment','%s_gene'%(EPI_target_label)]
	if args.conservative:
		use_cols = ['%s_gene'%(EPI_target_label)]
	my_query.get_final_target_genes(use_cols,target_gene_set_label)
	my_query.deg_evidence_filter(target_gene_set_label,args.deg_tsv,args.FDR_col_name,args.LFC_col_name,args.FDR_cutoff, args.LFC_cutoff)
	my_query.save_filtered_deg_table(target_gene_set_label)
	my_query.print_log()
	
	my_query.df.to_csv(args.output+".tsv",sep="\t")
	
	# python assign_targets.py -q example/NFIX_idr_peaks.bed -tss data/mm9.ensembl_v67.TSS.gene_name.bed -epi data/captureC.HSC.mm9.bed

	
if __name__ == "__main__":
	main()




































