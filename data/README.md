
## Preprocessing steps to get formatted captureC data

capture C data from: wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119339/suppl/GSE119339_SignificantInteractions_LiMACC.csv.gz


```

python get_mouse_captureC.py

python get_HPC7_captureC.py


There some genes not found in the TSS annotation file
>>> set(df.gene)-set(df1[3])
{'AI662270', 'RP23-3M10.11', 'Thoc4', 'E130006D01Rik', 'RP23-3M10.7', 'Xist', 'Neat1', 'Trfr2', 'Gas5', 'Tcfec', '4833418N02Rik', 'Igf2as', 'Bbip1', 'H19', '2810008D09Rik', '5830405N20Rik', 'Nlrc5', 'Tubb2c', '2310026L22Rik'}


```

## Files

TSS annotation: mm9.ensembl_v67.TSS.gene_name.bed

HSC captureC formatted: captureC.HSC.mm9.bed

HPC7 captureC formatted: HPC7.mm9.captureC.bed

mouse motif pwm: mouse_TF.meme


