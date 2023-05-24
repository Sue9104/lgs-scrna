#!/bin/sh
ANNOTATION=$1
MATRIX=$2
DESIGN=$3
OUTDIR=$4

TAPPAS_DIR=`dirname $0`/scripts/
# Prepare
[ ! -d ${OUTDIR} ] && mkdir ${OUTDIR}
[ ! -d ${OUTDIR}/InputData ] && mkdir ${OUTDIR}/InputData
cp ${ANNOTATION} ${OUTDIR}/InputData/annotations.gff3
sed -i '1i seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes' ${OUTDIR}/InputData/annotations.gff3
cp ${MATRIX} ${OUTDIR}/InputData/input_matrix.tsv
cp ${DESIGN} ${OUTDIR}/InputData/design.tsv
awk '$3=="transcript"' ${OUTDIR}/InputData/annotations.gff3|cut -f1,5> ${OUTDIR}/InputData/transcript_length.orig.tsv
Rscript ${TAPPAS_DIR}/design.R ${OUTDIR}/InputData/design.tsv ${OUTDIR}/Data/
Rscript ${TAPPAS_DIR}/PCA.R ${OUTDIR}/InputData/input_matrix.tsv ${OUTDIR}/Data/exp_factors.txt ${OUTDIR}/InputData/matrix_pca.tsv
Rscript ${TAPPAS_DIR}/tappas_inpMatrix.R ${OUTDIR}/InputData/input_matrix.tsv ${OUTDIR}/Data/exp_factors.txt ${OUTDIR}/InputData/transcript_length.orig.tsv ${OUTDIR}/InputData/input_normalized_matrix.tsv 1.0 100.0 Y

# Preprocess
[ ! -d ${OUTDIR}/Data ] && mkdir ${OUTDIR}/Data
Rscript ${TAPPAS_DIR}/preprocess.R ${OUTDIR}/InputData/annotations.gff3 ${OUTDIR}/InputData/design.tsv ${OUTDIR}/InputData/input_matrix.tsv ${OUTDIR}/InputData/input_normalized_matrix.tsv ${OUTDIR}/Data
Rscript ${TAPPAS_DIR}/PCA.R ${OUTDIR}/InputData/input_normalized_matrix.tsv ${OUTDIR}/Data/exp_factors.txt ${OUTDIR}/Data/matrix_pca.tsv

# DEA
Rscript ${TAPPAS_DIR}/DEA.R -d gene -m NOISEQ -r biological ${OUTDIR}/Data ${OUTDIR}/DEA
Rscript ${TAPPAS_DIR}/DEA.R -d trans -m NOISEQ -r biological ${OUTDIR}/Data ${OUTDIR}/DEA
Rscript ${TAPPAS_DIR}/DEA.R -d protein -m NOISEQ -r biological ${OUTDIR}/Data ${OUTDIR}/DEA

# DIU
Rscript ${TAPPAS_DIR}/DIU.R -d trans -m DEXSEQ -f 1.5 -t FOLD ${OUTDIR}/Data ${OUTDIR}/DIU
Rscript ${TAPPAS_DIR}/DIU_ScatterPlot.R -d trans ${OUTDIR}/DIU/trans_diu_matrix.tsv ${OUTDIR}/DIU/trans_scatter_plot.png
Rscript ${TAPPAS_DIR}/DIU.R -d protein -m DEXSEQ -f 1.5 -t FOLD ${OUTDIR}/Data ${OUTDIR}/DIU
Rscript ${TAPPAS_DIR}/DIU_ScatterPlot.R -d protein ${OUTDIR}/DIU/protein_diu_matrix.tsv ${OUTDIR}/DIU/protein_scatter_plot.png

# DFI
#Rscript ${TAPPAS_DIR}/DFI.R  -a908246893 -mDEXSEQ -i${OUTDIR}/Data -d${OUTDIR}/Data/DFI -o${OUTDIR}/DFI -f0 -tFOLD -s0.05 -c${OUTDIR}/Data/Content/dfi_total_features.908246893.tsv -g1${OUTDIR}/Data/Content/dfi_test_features.908246893.tsv -g2${OUTDIR}/Data/Content/dfi_test_genes.908246893.tsv -x${OUTDIR}/Data/Content/dfi_matrix.908246893.tsv -ltNMD,RBP_Binding,PAS,3UTRmotif,repeat,5UTR_Length,CDS_Length,polyA_Site,3UTR_Length,uORF,5UTRmotif,miRNA_Binding -lpSIGNAL,DOMAIN,TRANSMEM,ACT_SITE,INTRAMEM,TRANSMEM,PTM,BINDING,MOTIF,COILED,COMPBIAS,COILED,MOTIF,DISORDER
#Rscript ${TAPPAS_DIR}/DFI_BarPlot.R -i${OUTDIR}/Data/Content/dfi_matrix.908246893.tsv -a908246893 -c${OUTDIR}/Data/Content/dfi_total_features.908246893.tsv -ltNMD,RBP_Binding,PAS,3UTRmotif,repeat,5UTR_Length,CDS_Length,polyA_Site,3UTR_Length,uORF,5UTRmotif,miRNA_Binding -lpSIGNAL,DOMAIN,TRANSMEM,ACT_SITE,INTRAMEM,TRANSMEM,PTM,BINDING,MOTIF,COILED,COMPBIAS,COILED,MOTIF,DISORDER -o1${OUTDIR}/DFI/dfi_bar_plot1.908246893.png -o2${OUTDIR}/DFI/dfi_bar_plot2.908246893.png -o3${OUTDIR}/DFI/dfi_bar_plot3.908246893.png -o4${OUTDIR}/DFI/dfi_bar_plot4.908246893.png -t1${OUTDIR}/Data/Content/dfi_test_features.908246893.tsv -t2${OUTDIR}/Data/Content/dfi_test_genes.908246893.tsv

# DPA
Rscript ${TAPPAS_DIR}/DPA.R -m DEXSEQ -f 1.5 -t FOLD -l 60 ${OUTDIR}/InputData/design.tsv ${OUTDIR}/Data ${OUTDIR}/DPA
Rscript ${TAPPAS_DIR}/DPA_Heatmap.R -s 0.05 ${OUTDIR}/DPA ${OUTDIR}/DPA/heatmap_plot.png

# UTRL
Rscript ${TAPPAS_DIR}/UTRL.R -m WILCOXON -f 1.5 -t FOLD -l 100 ${OUTDIR}/InputData/design.tsv ${OUTDIR}/Data ${OUTDIR}/UTRL

#
#Rscript ${TAPPAS_DIR}/FEA.R -mWallenius -s2000 -u0 -dgene -f01756168879 -i${OUTDIR}/Data/FEA -o${OUTDIR}/FEA -g/home/minsu/tappasWorkspace/App/goAncestors.obo
#Rscript ${TAPPAS_DIR}/GSEA.R -dgene -mGOGLM -a236641335 -i${OUTDIR}/Data/GSEA -o${OUTDIR}/GSEA -vFALSE -g/home/minsu/tappasWorkspace/App/goAncestors.obo

# Stat
Rscript ${TAPPAS_DIR}/stat.R ${OUTDIR} ${OUTDIR}/Stat
