import pandas as pd
file1 = "mm9.ensembl_v67.gene_name.bed"
file2 = "captureC.HSC.mm9.bed\"
file2 = "captureC.HSC.mm9.bed"
df1 = pd.read_csv(file1,sep="\t",header=None)
df1.head()
df2 = pd.read_csv(file2,sep="\t",header=None)
df2.head()
set(df2[3])-set(df1[3])
import pandas as pd
file="GSE119339_SignificantInteractions_LiMACC.csv.gz"
df = pd.read_csv(file)
df.head()
df.columns
df.score = df[['score.BM','score.FL']].max()
df.head()
df[['score.BM','score.FL']].max().shape
df.score = df[['score.BM','score.FL']].max(axis=1)
df.gene = df.gene_name.apply(lambda x:x.split("_")[0])
file1 = "../example/mm9.ensembl_v67.gene_name.bed"
df1 = pd.read_csv(file1,sep="\t",header=None)
df1.head()
set(df.gene)-set(df[3])
set(df.gene)-set(df1[3])
df.head()
df.columns
df.chr = "chr"+df.otherEnd_chr
df[['chr','otherEnd_start','otherEnd_end','name','score']].to_csv("captureC.HSC.mm9.bed",sep="\t",header=False,index=False)
df.gene = df.gene_name.apply(lambda x:x.split("_")[0])
df.columns
df.gene.head()
df['gene'] = df.gene_name.apply(lambda x:x.split("_")[0])
df['score'] = df[['score.BM','score.FL']].max(axis=1)
df['chr'] = "chr"+df.otherEnd_chr
df[['chr','otherEnd_start','otherEnd_end','name','score']].to_csv("captureC.HSC.mm9.bed",sep="\t",header=False,index=False)
df[['chr','otherEnd_start','otherEnd_end','gene','score']].to_csv("captureC.HSC.mm9.bed",sep="\t",header=False,index=False)
