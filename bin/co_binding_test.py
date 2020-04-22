#!/usr/bin/env python


"""

Given two TFBS bed file, check overlap and output p-value

Parameters
-----------

d distance between the two peaks to consider as co binding

Dependency
-----------

bedtools

Output
------

Venn, FG overlap, mean BG overlap, p-value

"""


from tf_target_finder.utils import *
from scipy.stats import chi2_contingency


def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	mainParser.add_argument('-f1',"--query",help="Query bed. 3 column bed file, additional columns are OK, but will be ignored",required=True)	
	mainParser.add_argument('-f2',"--co_tf",help="3 column bed file, additional columns are OK, but will be ignored",required=True)	
	mainParser.add_argument('-bg',help="query TF all peaks as BG",required=True)	
	
	mainParser.add_argument("-d",help="distance cutoff for overlap",default=0,type=int)	


	mainParser.add_argument('-o',"--output",  help="output name",default="co_binding_test")
	mainParser.add_argument('-g',"--chrom_size",  help="output name",default="/home/yli11/Data/Mouse/mm9/annotations/mm9_main.chrom.sizes")
	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args

def sort_bed(input,output):
	os.system("sort -k1,1 -k2,2n %s>%s"%(input,output))

def bedtools_closest(f1,f2,d):
	outFile = str(uuid.uuid4()).split("-")[-1]
	sort1 = str(uuid.uuid4()).split("-")[-1]
	sort2 = str(uuid.uuid4()).split("-")[-1]
	sort_bed(f1,sort1)
	sort_bed(f2,sort2)
	command = "bedtools closest -a %s -b %s -d > %s"%(sort1,sort2,outFile)
	os.system(command)
	df = pd.read_csv(outFile,sep="\t",header=None)
	# print (df.head())
	# print (df.shape)
	delete_file(outFile)
	delete_file(sort1)
	delete_file(sort2)
	df = df[df[df.columns[-1]]<=d]
	# print (df.head())
	if df.shape[0] == 0:
		return 0
	df['name'] = df[0]+":"+df[1].astype(str)+"-"+df[2].astype(str)
	df = df.drop_duplicates('name')
	return df.shape[0]
	
def get_background_list(query_bed,chrom_size,co_df_bed,d):
	outFile = str(uuid.uuid4()).split("-")[-1]
	bedtools_shuffle = "bedtools shuffle -seed 1 -i %s -excl %s -g %s -chrom > %s"%(query_bed,query_bed,chrom_size,outFile)
	os.system(bedtools_shuffle)
	count = bedtools_closest(outFile,co_df_bed,d)
	delete_file(outFile)
	return count
	
def z_test_p_value(fg,bg_list):
	mean = np.mean(bg_list)
	std = np.std(bg_list)
	z_score = (fg-mean)/std
	return mean,scipy.stats.norm.sf(abs(z_score))*2

def chi2_test(FG_count,FG_overlap,BG_count,BG_overlap):
	A = FG_overlap
	B = FG_count - FG_overlap
	C = BG_overlap
	D = BG_count - BG_overlap
	odd,p_value,t1,t2 = chi2_contingency([[A,B],[C,D]])
	return p_value

def main():

	args = my_args()
	query = read_bed(args.query,0)
	tf = read_bed(args.co_tf,0)
	query['name'] = query[0]+":"+query[1].astype(str)+"-"+query[2].astype(str)
	query = query.drop_duplicates('name')
	FG_count = query.shape[0]
	print ("Input number of peaks in query bed: %s"%(FG_count))
	print ("----------- %s --------------"%(args.co_tf))
	FG_overlap = bedtools_closest(args.query,args.co_tf,args.d)
	FG_percent = FG_overlap/float(FG_count)
	print ("Number of TFBSs overlaped with query: {0} ({1:.1%})".format(FG_overlap,FG_percent))
	
	
	bg_query = read_bed(args.bg,0)
	bg_query['name'] = bg_query[0]+":"+bg_query[1].astype(str)+"-"+bg_query[2].astype(str)
	bg_query = bg_query.drop_duplicates('name')
	BG_count = bg_query.shape[0] - FG_count
	print ("Input number of peaks in background bed: %s"%(BG_count))
	BG_overlap = bedtools_closest(args.bg,args.co_tf,args.d) - FG_overlap
	print ("Number of TFBSs overlaped with background: %s"%(BG_overlap))
	
	P_value = chi2_test(FG_count,FG_overlap,BG_count,BG_overlap)
	
	print ("Overlapping p_value: %s. Enrichment: %s."%(P_value,FG_percent/(float(BG_overlap)/BG_count)))
	

	

	
	
if __name__ == "__main__":
	main()

































