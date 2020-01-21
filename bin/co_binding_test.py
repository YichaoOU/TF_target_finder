#!/usr/bin/env python


"""

Given two TFBS bed file, check overlap and output p-value

Parameters
-----------

d1   extend of f1
d2   extend of f2

Dependency
-----------

bedtools

"""


import sys
import os
p_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
sys.path.append(os.path.abspath(p_dir+"../utils/"))
from utils import *



def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	mainParser.add_argument('-f1',"--bed1",help="Query bed. 3 column bed file, additional columns are OK, but will be ignored",required=True)	
	mainParser.add_argument('-f2',"--bed2",help="3 column bed file, additional columns are OK, but will be ignored",required=True)	
	
	mainParser.add_argument("-d1",help="distance cutoff for TF bed",default=0,type=int)	
	mainParser.add_argument("-d2",help="distance cutoff for Gene-associated features bed file",default=0,type=int)	

	mainParser.add_argument('-o',"--output",  help="output overlapped gene list and peak list, [x].gene.list and [x].TFBS.list, to find out which peak overlaps which gene, see gene_TF.intersect.[x].bed.",default="output")
	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args



	
def main():

	args = my_args()
	tf1 = read_bed(args.bed1,args.d1)
	tf2 = read_bed(args.bed2,args.d2)

	
	
if __name__ == "__main__":
	main()

































