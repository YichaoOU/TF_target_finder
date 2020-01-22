

################### PARAMETERS ##########################
genome="mm9"
logFC=1
FDR=0.05
TSS_flank=10000
peak_flank=0
capC_flank=0
output_label="NFIX"

co_binding_distance=500

co_motif_distance=200

motif_scanning_flank=200

################### INPUTS ##############################
query="NFIX_idr_peaks.bed"
EPI="captureC.HSC.mm9.bed"
TSS="mm9.ensembl_v67.gene_name.bed"
gene_expression="results.KO_vs_WT.txt"
BIOGRID_mouse_known_interaction="Mus_musculus.interaction.list"
BIOGRID_ALL_known_interaction="ALL.interaction.list"


################### MAIN ################################

echo "Extracting targeted genes based on enhancers"

# 1. TF overlap enhancer

../bin/TF_overlap.py -f1 $query -f2 $EPI -d1 $peak_flank -d2 $capC_flank -o $output_label.capC

echo "Extracting targeted genes based on promoters"

# 2. TF overlap TSS

../bin/TF_overlap.py -f1 $query -f2 $TSS -d1 $peak_flank -d2 $TSS_flank -o $output_label.TSS

# Union of step 1 and step 2

echo "Combining gene list from promoters and enhancers"

cat $output_label.capC.gene.list $output_label.TSS.gene.list > $output_label.gene.list

cat $output_label.capC.TFBS.list $output_label.TSS.TFBS.list > $output_label.TFBS.list


# 3. gene expression evidence

echo "Filter targeted genes based on gene expression from KO experiments"

../bin/evidence_filter.py -g $output_label.gene.list -e $gene_expression --cols gene,logFC,adj.P.Val --cutoff $logFC,$FDR -o $output_label.DEG.list --TFBS_list $output_label.TFBS.list

echo "Filter targeted genes based on known interactions"

../bin/evidence_filter.py -g $output_label.gene.list -e $BIOGRID_ALL_known_interaction -o $output_label.BIOGRID.ALL.list -t $output_label --TFBS_list $output_label.TFBS.list

echo "Output genes that are supported by either evidence"

cat $output_label.BIOGRID.ALL.gene.list $output_label.DEG.gene.list > $output_label.direct_targets.gene.list

cat $output_label.BIOGRID.ALL.TFBS.list $output_label.DEG.TFBS.list > $output_label.direct_targets.TFBS.list

echo "Final direct targets list is: $output_label.direct_targets.list"

# 4. co-binding test

echo "Begin co-binding test for chip-seq peaks"

cd co_factors_chip_seq

../../bin/co_binding_test.py -f1 ../$output_label.direct_targets.TFBS.list -f2 GSM1708650_Erg_416B.mm9.bed -d $co_binding_distance -bg ../$query


../../bin/co_binding_test.py -f1 ../$output_label.direct_targets.TFBS.list -f2 GSM1708655_Meis1_416B.mm9.bed -d $co_binding_distance  -bg ../$query

../../bin/co_binding_test.py -f1 ../$output_label.direct_targets.TFBS.list -f2 GSM1708651_Fli1_416B.mm9.bed -d $co_binding_distance  -bg ../$query

../../bin/co_binding_test.py -f1 ../$output_label.direct_targets.TFBS.list -f2 GSM1708656_Sfpi1_416B.mm9.bed -d $co_binding_distance  -bg ../$query

../../bin/co_binding_test.py -f1 ../$output_label.direct_targets.TFBS.list -f2 GSM1708652_Gata2_416B.mm9.bed -d $co_binding_distance  -bg ../$query

../../bin/co_binding_test.py -f1 ../$output_label.direct_targets.TFBS.list -f2 GSM1708657_Runx1_416B.mm9.bed -d $co_binding_distance  -bg ../$query

../../bin/co_binding_test.py -f1 ../$output_label.direct_targets.TFBS.list -f2 GSM1708653_Gfi1b_416B.mm9.bed -d $co_binding_distance  -bg ../$query

../../bin/co_binding_test.py -f1 ../$output_label.direct_targets.TFBS.list -f2 GSM1708658_Scl_416B.mm9.bed -d $co_binding_distance  -bg ../$query

../../bin/co_binding_test.py -f1 ../$output_label.direct_targets.TFBS.list -f2 GSM1708654_Lyl1_416B.mm9.bed -d $co_binding_distance  -bg ../$query

# co-binding motif test

cd ..

cd co_factors_motif

echo "Begin co-binding test for known motifs"

echo "Perform motif scanning"


../../bin/motif_scanning.py -f ../$query -m NFIX_mouse_known_motifs.meme -o NFIX_motif.BG.match -e $motif_scanning_flank

../../bin/motif_scanning.py -f ../$output_label.direct_targets.TFBS.list -m NFIX_mouse_known_motifs.meme -o NFIX_motif.FG.match -e $motif_scanning_flank

../../bin/motif_scanning.py -f ../$query -m mouse_TF.meme -o ERG_motif.match -e $motif_scanning_flank --motif_ids Erg_MA0474.1,ERG_MOUSE.H11MO.0.A

../../bin/motif_scanning.py -f ../$query -m mouse_TF.meme -o GATA2_motif.match -e $motif_scanning_flank --motif_ids GATA2_MOUSE.H11MO.0.A

../../bin/motif_scanning.py -f ../$query -m mouse_TF.meme -o RUNX1_motif.match -e $motif_scanning_flank --motif_ids RUNX1_MA0002.2,RUNX1_MOUSE.H11MO.0.A

../../bin/motif_scanning.py -f ../$query -m mouse_TF.meme -o GFI1B_motif.match -e $motif_scanning_flank --motif_ids Gfi1b_MA0483.1,GFI1B_MOUSE.H11MO.0.A

../../bin/motif_scanning.py -f ../$query -m mouse_TF.meme -o Sfpi1_motif.match -e $motif_scanning_flank --motif_ids Sfpi1_M6122_1.02,Sfpi1_primary_UP00085_1,Sfpi1_secondary_UP00085_2

../../bin/motif_scanning.py -f ../$query -m mouse_TF.meme -o MEIS1_motif.match -e $motif_scanning_flank --motif_ids MEIS1_MOUSE.H11MO.0.A,Meis1_MA0498.1

../../bin/motif_scanning.py -f ../$query -m mouse_TF.meme -o TAL1_motif.match -e $motif_scanning_flank --motif_ids TAL1_MOUSE.H11MO.0.A,Tal1_Gata1_MA0140.1

../../bin/motif_scanning.py -f ../$query -m mouse_TF.meme -o LYL1_motif.match -e $motif_scanning_flank --motif_ids LYL1_MOUSE.H11MO.0.A

../../bin/motif_scanning.py -f ../$query -m mouse_TF.meme -o FLI1_motif.match -e $motif_scanning_flank --motif_ids FLI1_MOUSE.H11MO.0.A


../../bin/co_binding_test.py -f1 NFIX_motif.FG.match.bed -f2 ERG_motif.match.bed -d $co_motif_distance -bg NFIX_motif.BG.match.bed
../../bin/co_binding_test.py -f1 NFIX_motif.FG.match.bed -f2 GATA2_motif.match.bed -d $co_motif_distance -bg NFIX_motif.BG.match.bed
../../bin/co_binding_test.py -f1 NFIX_motif.FG.match.bed -f2 RUNX1_motif.match.bed -d $co_motif_distance -bg NFIX_motif.BG.match.bed
../../bin/co_binding_test.py -f1 NFIX_motif.FG.match.bed -f2 GFI1B_motif.match.bed -d $co_motif_distance -bg NFIX_motif.BG.match.bed
../../bin/co_binding_test.py -f1 NFIX_motif.FG.match.bed -f2 Sfpi1_motif.match.bed -d $co_motif_distance -bg NFIX_motif.BG.match.bed
../../bin/co_binding_test.py -f1 NFIX_motif.FG.match.bed -f2 MEIS1_motif.match.bed -d $co_motif_distance -bg NFIX_motif.BG.match.bed
../../bin/co_binding_test.py -f1 NFIX_motif.FG.match.bed -f2 TAL1_motif.match.bed -d $co_motif_distance -bg NFIX_motif.BG.match.bed
../../bin/co_binding_test.py -f1 NFIX_motif.FG.match.bed -f2 LYL1_motif.match.bed -d $co_motif_distance -bg NFIX_motif.BG.match.bed
../../bin/co_binding_test.py -f1 NFIX_motif.FG.match.bed -f2 FLI1_motif.match.bed -d $co_motif_distance -bg NFIX_motif.BG.match.bed


