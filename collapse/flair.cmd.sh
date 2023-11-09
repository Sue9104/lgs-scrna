basedir=$(dirname "$0")
inbam=$1
infa=$2
intsv=$3
outdir=$4
prefix=$5
genome=$6
gtf=$7
#genome=/public/home/msu/genomes/hg38/ensembl/Homo_sapiens.GRCh38.dna_sm.primary_assembly.chr.fa
#gtf=/public/home/msu/genomes/hg38/gencode/gencode.v42.primary_assembly.annotation.gtf
cage=/public/home/msu/pipelines/exon-usages-in-3rd-seq/SQANTI3/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed

cores=20
# create outdir if not exists
[ ! -d "${outdir}" ] && mkdir ${outdir}

# bam to bed12
echo "samtools index ${inbam} && bam2Bed12 -i ${inbam} > ${outdir}/${prefix}.bed"
samtools index ${inbam}
[ ! -f "${outdir}/${prefix}.bed" ] && bam2Bed12 -i ${inbam} > ${outdir}/${prefix}.bed
cut -f1-3 ${cage} > ${outdir}/promoter_regions.bed

# correct
echo flair correct -q ${outdir}/${prefix}.bed -f ${gtf} -g ${genome} -o ${outdir}/${prefix}
[ ! -f "${outdir}/${prefix}.corrected_all_corrected.bed" ] && flair correct -q ${outdir}/${prefix}.bed -f ${gtf} -g ${genome} -o ${outdir}/${prefix}.corrected

# collapse
echo flair collapse -r ${infa} -q ${outdir}/${prefix}.corrected_all_corrected.bed -g ${genome} -t ${cores} --generate_map --temp_dir ${outdir}/temp_collapse --keep_intermediate -f ${gtf} -o ${outdir}/${prefix}.collapse -p ${outdir}/promoter_regions.bed --trust_ends --no_gtf_end_adjustment --isoformtss --max_ends 5
[ ! -f "${outdir}/${prefix}.collapse.isoforms.gtf" ] && flair collapse -r ${infa} -q ${outdir}/${prefix}.corrected_all_corrected.bed -g ${genome} -t ${cores} --generate_map --temp_dir ${outdir}/temp_collapse --keep_intermediate -f ${gtf} -o ${outdir}/${prefix}.collapse -p ${outdir}/promoter_regions.bed --trust_ends --no_gtf_end_adjustment --isoformtss --max_ends 5

#	quantify
echo flair quantify -r ${intsv} -i ${outdir}/${prefix}.collapse.isoforms.fa --isoform_bed ${outdir}/${prefix}.collapse.isoforms.bed --generate_map --tpm --check_splice --stringent -o ${outdir}/${prefix}.quantify
flair quantify -r ${intsv} -i ${outdir}/${prefix}.collapse.isoforms.fa --isoform_bed ${outdir}/${prefix}.collapse.isoforms.bed --generate_map --tpm --check_splice --stringent -o ${outdir}/${prefix}.quantify

# diffexp
echo flair diffExp -q ${outdir}/${prefix}.quantify.counts.tsv -o ${outdir}/diffexp -e 1 -of
#flair diffExp -q ${outdir}/${prefix}.quantify.counts.tsv -o ${outdir}/diffexp -e 1 -of

# diffsplice
echo flair diffSplice -i ${outdir}/${prefix}.collapse.isoforms.bed -q ${outdir}/${prefix}.quantify.counts.tsv --test -o ${outdir}/diffsplice -of
#flair diffSplice -i ${outdir}/${prefix}.collapse.isoforms.bed -q ${outdir}/${prefix}.quantify.counts.tsv --test -o ${outdir}/diffsplice -of

