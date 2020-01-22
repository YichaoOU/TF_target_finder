

################### PARAMETERS ##########################
genome="mm9"
logFC=2
FDR=0.05
TSS_flank=10000
peak_flank=0
capC_flank=0
output_label="NFIX"



################### INPUTS ##############################
query="NFIX_idr_peaks.bed"
EPI="captureC.HSC.mm9.bed"
TSS="mm9.ensembl_v67.gene_name.bed"
gene_expression="results.KO_vs_WT.txt"
BIOGRID_mouse_known_interaction="Mus_musculus.interaction.list"
BIOGRID_ALL_known_interaction="ALL.interaction.list"


################### MAIN ################################

# 1. TF overlap enhancer

../bin/TF_overlap.py -f1 $query -f2 $EPI -d1 $peak_flank -d2 $capC_flank -o $output_label.capC

# 2. TF overlap TSS

../bin/TF_overlap.py -f1 $query -f2 $TSS -d1 $peak_flank -d2 $TSS_flank -o $output_label.TSS

# Union of step 1 and step 2

cat $output_label.capC.gene.list $output_label.TSS.gene.list > $output_label.gene.list

cat $output_label.capC.TFBS.list $output_label.TSS.TFBS.list > $output_label.TFBS.list


# 3. gene expression evidence

../bin/evidence_filter.py -g $output_label.gene.list -e $gene_expression --cols gene,logFC,adj.P.Val --cutoff $logFC,$FDR -o $output_label.DEG.list

../bin/evidence_filter.py -g $output_label.gene.list -e $BIOGRID_ALL_known_interaction -o $output_label.BIOGRID.ALL.list -t $output_label


../../bin/co_binding_test.py -f1 ../NFIX.TFBS.list -f2 GSM1708650_Erg_416B.mm9.bed -d 1000

