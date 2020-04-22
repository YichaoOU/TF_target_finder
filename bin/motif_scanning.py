#!/usr/bin/env python


"""


Dependency
----------

Bedtools
MEME

Input
-----

chip-seq peak file
motif meme file

Output
-----

motif occurrence bed file

Parameters
---------

-extend flank length of the peak file

"""


from tf_target_finder.utils import *
from utils import *



def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	mainParser.add_argument('-f',"--TF_bed",help="3 column bed file, additional columns are OK, but will be ignored",required=True)	
	mainParser.add_argument('-m',"--motif_meme",help="motif meme file",required=True)	
	mainParser.add_argument('--motif_ids',help="input motif ids sep by ,",default="None")	
	
	mainParser.add_argument("-e",'--extend',help="extend search on the flank sequences",default=0,type=int)	
	mainParser.add_argument("-fa",'--genome_fasta',help="genome fasta sequence",default="/home/yli11/Data/Mouse/mm9/index/STAR/Mus_musculus.NCBIM37.67.add_chr.dna.toplevel.fa")	


	mainParser.add_argument('-o',"--output",  help="output file name",default="output")
	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args

def main():

	args = my_args()
	## extend bed
	tf = read_bed(args.TF_bed,args.extend)
	tmp_file = str(uuid.uuid4()).split("-")[-1]
	tmp_fa = str(uuid.uuid4()).split("-")[-1]
	to_bed(tf[[0,1,2]],tmp_file)
	
	## get fasta
	command = "bedtools getfasta -fi %s -fo %s -bed %s"%(args.genome_fasta,tmp_fa,tmp_file)
	os.system(command)
	delete_file(tmp_file)
	
	## motif scan
	motif_ids = ""
	command = "fimo {{motif_ids}} --text --verbosity 1 --parse-genomic-coord %s %s > {{output}};sed -i '1d' {{output}};cut -f 2,3,4 {{output}} > %s.bed"%(args.motif_meme,tmp_fa,args.output)
	tmp_fimo = str(uuid.uuid4()).split("-")[-1]
	command = command.replace("{{output}}",tmp_fimo)
	if args.motif_ids != "None":
		for i in args.motif_ids.split(","):
			motif_ids+="--motif %s "%(i)
	command = command.replace("{{motif_ids}}",motif_ids)
	os.system(command)
	delete_file(tmp_fa)
	delete_file(tmp_fimo)
	
	
if __name__ == "__main__":
	main()

































