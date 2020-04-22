import pandas as pd
df = pd.read_csv("HPC7_Promoter_Capture_Interactions.ibed",sep="\t")

myList = [x.split(", ") for x in df.bait_name]
myList = [item for sublist in myList for item in sublist]
myList = list(set(myList))

my2 = [x.split("-")[-1] for x in myList]
my2=list(set(my2))
my2.remove("Car1")
my2.remove("uce")

def define_gene_name(x):
	genes = x.split(",")
	new_genes = []
	special_char = ['Car1','uce']
	for g in genes:
		if g.split("-")[-1] in special_char:
			new_genes.append(g)
		else:
			new_genes.append("-".join(g.split("-")[:-1]))
	return list(set(new_genes))
	
df['gene_name'] = [define_gene_name(x) for x in df.bait_name]

df = df.explode('gene_name')
cols=['otherEnd_chr','otherEnd_start','otherEnd_end','gene_name','score']
df['otherEnd_chr'] = 'chr'+df['otherEnd_chr']
df[cols].to_csv("HPC7.mm9.captureC.bed",index=False,header=False,sep="\t")

file1 = "../example/mm9.ensembl_v67.gene_name.bed"
df1 = pd.read_csv(file1,sep="\t",header=None)

print (set(df.gene_name)-set(df1[3]))


 # otherEnd_chr  otherEnd_start score
 # otherEnd_chr  otherEnd_start score