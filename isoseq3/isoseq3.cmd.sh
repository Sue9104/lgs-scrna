basedir=$(dirname "$0")
input=$1
outdir=$2
prefix=$3

genome=/public/home/msu/genomes/hg38/ensembl/Homo_sapiens.GRCh38.dna_sm.primary_assembly.chr.fa
gtf=/public/home/msu/genomes/hg38/gencode/gencode.v42.primary_assembly.annotation.gtf
cage=/public/home/msu/pipelines/exon-usages-in-3rd-seq/SQANTI3/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed
polya=/public/home/msu/pipelines/exon-usages-in-3rd-seq/SQANTI3/data/polyA_motifs/mouse_and_human.polyA_motif.txt
primers=/public/home/msu/projects/seq3/examples/primers.fasta
barcodes=/public/home/msu/projects/seq3/examples/3M-february-2018-REVERSE-COMPLEMENTED.txt.gz
python3=/public/home/msu/miniconda3/bin/python3
# create outdir if not exists
[ ! -d "${outdir}" ] && mkdir ${outdir}

# cDNA primer removal and read orientation
echo lima --per-read --isoseq ${input} ${primers} ${outdir}/${prefix}.output.bam
[ ! -f "${outdir}/${prefix}.output.5p--3p.bam" ] && lima --per-read --isoseq ${input} ${primers} ${outdir}/${prefix}.output.bam

# Clip UMI and cell barcode
echo isoseq3 tag ${outdir}/${prefix}.output.5p--3p.bam ${outdir}/${prefix}.flt.bam --design T-12U-16B
[ ! -f "${outdir}/${prefix}.flt.bam" ] && isoseq3 tag ${outdir}/${prefix}.output.5p--3p.bam ${outdir}/${prefix}.flt.bam --design T-12U-16B

# Remove poly(A) tails and concatemer
echo isoseq3 refine ${outdir}/${prefix}.flt.bam ${primers} ${outdir}/${prefix}.fltnc.bam --require-polya
[ ! -f "${outdir}/${prefix}.fltnc.bam" ] && isoseq3 refine ${outdir}/${prefix}.flt.bam ${primers} ${outdir}/${prefix}.fltnc.bam --require-polya

# Correct single cell barcodes based on an include list
echo isoseq3 correct -B ${barcodes} ${outdir}/${prefix}.fltnc.bam ${outdir}/${prefix}.corrected.bam
[ ! -f "${outdir}/${prefix}.corrected.bam" ] && isoseq3 correct -B ${barcodes} ${outdir}/${prefix}.fltnc.bam ${outdir}/${prefix}.corrected.bam
# Barcode Statistics Documentation
echo isoseq3 bcstats --json ${outdir}/${prefix}.bcstats.json -o ${outdir}/${prefix}.bcstats.tsv ${outdir}/${prefix}.corrected.bam
[ ! -f "${outdir}/${prefix}.bcstats.tsv" ] && isoseq3 bcstats --json ${outdir}/${prefix}.bcstats.json -o ${outdir}/${prefix}.bcstats.tsv ${outdir}/${prefix}.corrected.bam
echo ${python3} ${basedir}/plot_knees.py -t ${outdir}/${prefix}.bcstats.tsv -o ${outdir}/${prefix}.cellbackground --estimate_percentile 95
[ ! -f "${outdir}/${prefix}.cellbackground.knee.png" ] && ${python3} ${basedir}/plot_knees.py -t ${outdir}/${prefix}.bcstats.tsv -o ${outdir}/${prefix}.cellbackground --estimate_percentile 95

# Deduplicate reads based on UMIs
echo samtools sort -@ 40 -t CB ${outdir}/${prefix}.corrected.bam -o ${outdir}/${prefix}.corrected.sorted.bam
[ ! -f "${outdir}/${prefix}.corrected.sorted.bam" ] && samtools sort -@ 40 -t CB ${outdir}/${prefix}.corrected.bam -o ${outdir}/${prefix}.corrected.sorted.bam
echo isoseq3 groupdedup ${outdir}/${prefix}.corrected.sorted.bam ${outdir}/${prefix}.dedup.bam
[ ! -f "${outdir}/${prefix}.dedup.bam" ] && isoseq3 groupdedup ${outdir}/${prefix}.corrected.sorted.bam ${outdir}/${prefix}.dedup.bam

# Map reads to a reference genom
echo pbmm2 align --preset ISOSEQ --sort ${outdir}/${prefix}.dedup.bam ${genome} ${outdir}/${prefix}.aligned.bam
[ ! -f "${outdir}/${prefix}.aligned.bam" ] && pbmm2 align --preset ISOSEQ --sort ${outdir}/${prefix}.dedup.bam ${genome} ${outdir}/${prefix}.aligned.bam

# Collapse into unique isoforms
## Single-cell IsoSeq
echo isoseq3 collapse --keep-non-real-cells ${outdir}/${prefix}.aligned.bam ${outdir}/${prefix}.collapse.gff
[ ! -f "${outdir}/${prefix}.collapse.gff" ] && isoseq3 collapse --keep-non-real-cells ${outdir}/${prefix}.aligned.bam ${outdir}/${prefix}.collapse.gff
## Bulk IsoSeq
#isoseq3 collapse --do-not-collapse-extra-5exons mapped.bam collapsed.gff

# Sort input transcript GFF
echo pigeon sort ${outdir}/${prefix}.collapse.gff -o ${outdir}/${prefix}.sorted.gff
[ ! -f "${outdir}/${prefix}.sorted.gff" ] && pigeon sort ${outdir}/${prefix}.collapse.gff -o ${outdir}/${prefix}.sorted.gff

# Index the reference files
#pigeon index gencode.annotation.gtf
#pigeon index cage.bed
#pigeon index intropolis.tsv

# Classify Isoforms
#pigeon classify sorted.gff annotations.gtf reference.fa
echo pigeon classify ${outdir}/${prefix}.sorted.gff ${gtf} ${genome} --fl ${outdir}/${prefix}.collapse.abundance.txt --cage-peak ${cage} --poly-a ${polya} -d ${outdir} -o ${prefix}
[ ! -f "${outdir}/${prefix}_classification.txt" ] && pigeon classify ${outdir}/${prefix}.sorted.gff ${gtf} ${genome} --fl ${outdir}/${prefix}.collapse.abundance.txt --cage-peak ${cage} --poly-a ${polya} -d ${outdir} -o ${prefix}

# Filter isoforms
echo pigeon filter ${outdir}/${prefix}_classification.txt --isoforms ${outdir}/${prefix}.sorted.gff
[ ! -f "${outdir}/${prefix}_classification.filtered_lite_classification.txt" ] && pigeon filter ${outdir}/${prefix}_classification.txt --isoforms ${outdir}/${prefix}.sorted.gff

# Report gene saturation
echo pigeon report ${outdir}/${prefix}_classification.filtered_lite_classification.txt ${outdir}/${prefix}.saturation.txt
[ ! -f "${outdir}/${prefix}.saturation.txt" ] && pigeon report ${outdir}/${prefix}_classification.filtered_lite_classification.txt ${outdir}/${prefix}.saturation.txt


# Make Seurat compatible input
echo pigeon make-seurat --keep-ribo-mito-genes --dedup ${outdir}/${prefix}.dedup.fasta --group ${outdir}/${prefix}.collapse.group.txt -d ${outdir} -o ${prefix} ${outdir}/${prefix}_classification.filtered_lite_classification.txt
pigeon make-seurat --keep-ribo-mito-genes --dedup ${outdir}/${prefix}.dedup.fasta --group ${outdir}/${prefix}.collapse.group.txt -d ${outdir} -o ${prefix} ${outdir}/${prefix}_classification.filtered_lite_classification.txt
