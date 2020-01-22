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


import sys
import os
p_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
sys.path.append(os.path.abspath(p_dir+"../utils/"))
from utils import *



def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	mainParser.add_argument('-f1',"--TF_bed",help="3 column bed file, additional columns are OK, but will be ignored",required=True)	
	mainParser.add_argument('-f2',"--gene_bed",help="4 column bed file, the 4th column should be gene name. Additional columns are OK, but will be ignored",required=True)	
	
	mainParser.add_argument("-d1",help="distance cutoff for TF bed",default=0,type=int)	
	mainParser.add_argument("-d2",help="distance cutoff for Gene-associated features bed file",default=0,type=int)	

	mainParser.add_argument('-o',"--output",  help="output overlapped gene list and peak list, [x].gene.list and [x].TFBS.list, to find out which peak overlaps which gene, see gene_TF.intersect.[x].bed.",default="output")
	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args


def bedtools_overlap(tf,gene,output):
	f1 = str(uuid.uuid4()).split("-")[-1]
	f2 = str(uuid.uuid4()).split("-")[-1]
	to_bed(tf[[0,1,2]],f2)
	to_bed(gene[[0,1,2,3]],f1)
	command = "bedtools intersect -a %s -b %s -wo > gene_TF.intersect.%s.bed"%(f1,f2,output)
	os.system(command)
	delete_file(f1)
	delete_file(f2)
	return "gene_TF.intersect.%s.bed"%(output)
	
def main():

	args = my_args()
	tf = read_bed(args.TF_bed,args.d1)
	gene = read_bed(args.gene_bed,args.d2)
	overlaps = bedtools_overlap(tf,gene,args.output)
	print ("TF gene overlap file: %s"%(overlaps))
	df = read_bed(overlaps,0)
	tmp = df[[4,5,6,3]]
	tmp['name'] = tmp[4]+tmp[5].astype(str)+tmp[6].astype(str)+tmp[3].astype(str)
	tmp = tmp.drop_duplicates('name')
	tmp = tmp.drop(['name'],axis=1)
	to_bed(tmp,"%s.TFBS.list"%(args.output))
	print ("Overlapped TFBS list: %s.TFBS.list"%(args.output))
	df = df.drop_duplicates(3)
	to_bed(df[3],"%s.gene.list"%(args.output))
	print ("Overlapped gene list: %s.gene.list"%(args.output))
	
	
if __name__ == "__main__":
	main()

































