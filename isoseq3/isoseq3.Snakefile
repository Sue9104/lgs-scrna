import snakemake.utils
snakemake.utils.min_version('6.5.3')

DEFAULT_MEM_MB=64 * 1024  # 64 GB
DEFAULT_TIME_HOURS=24

import os
import pandas as pd
from dotenv import load_dotenv
load_dotenv()
mailto = os.getenv('mailto')
basedir = os.getenv('basedir')
conda_wrapper = os.path.join(basedir, "isoseq3.conda_wrapper")
print(conda_wrapper)

project=config['project']
genome=config['genome']
gtf=config['gtf']
info = pd.read_csv(config["infile"])
samples = info["sample"].tolist()

outdir=os.path.abspath(config['outdir'])
if not os.path.exists(outdir): os.makedirs(outdir)
for sample in samples:
    sample_dir = os.path.join(outdir, sample)
    if not os.path.exists(sample_dir): os.makedirs(sample_dir)


onsuccess:
    print('isoseq3 workflow success')
    #shell("echo Succeed~| mail -s 'Mito Calling Finished: {project}' {mailto}")
onerror:
    print('isoseq3 workflow error')
    #shell("echo Succeed~| mail -s 'Mito Calling Finished: {project}' {mailto}")


localrules: all
rule all:
    input:
        expand("{outdir}/{sample}/genes_seurat/matrix.mtx", outdir=outdir, sample = samples),

def get_bam_from_wildcards(wildcards):
    bam = info[info["sample"]== wildcards.sample]["bam"].squeeze(),
    return bam

rule remove_primers:
    output:
        "{outdir}/{sample}/{sample}.remove_primers.5p--3p.bam"
    log:
        out='{outdir}/{sample}/1_remove_primers_log.out',
        err='{outdir}/{sample}/1_remove_primers_log.err',
    params:
        inbam = get_bam_from_wildcards,
        primers = config['primers']
    resources:
        mem_mb=16 * 1024,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{conda_wrapper} lima --per-read --isoseq '
        ' {params.inbam}'
        ' {params.primers}'
        ' {outdir}/{sample}/{sample}.remove_primers.bam'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule tag:
    input:
        rules.remove_primers.output,
    output:
        "{outdir}/{sample}/{sample}.flt.bam",
    log:
        out='{outdir}/{sample}/2_tag_log.out',
        err='{outdir}/{sample}/2_tag_log.err',
    threads: 20
    resources:
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{conda_wrapper} isoseq3 tag'
        ' {input}'
        ' {output}'
        ' -j {threads}'
        ' --design T-12U-16B'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule refine:
    input:
        rules.tag.output
    output:
        "{outdir}/{sample}/{sample}.fltnc.bam",
    log:
        out='{outdir}/{sample}/3_refine_log.out',
        err='{outdir}/{sample}/3_refine_log.err',
    threads: 20
    params:
        primers = config['primers']
    resources:
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{conda_wrapper} isoseq3 refine'
        ' {input}'
        ' {params.primers}'
        ' {output}'
        ' -j {threads}'
        ' --require-polya'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule correct:
    input:
        rules.refine.output
    output:
        "{outdir}/{sample}/{sample}.correct.bam",
    log:
        out='{outdir}/{sample}/4_correct_log.out',
        err='{outdir}/{sample}/4_correct_log.err',
    threads: 20
    params:
        barcodes = config['barcodes']
    resources:
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{conda_wrapper} isoseq3 correct'
        ' -j {threads}'
        ' -B {params.barcodes}'
        ' {input}'
        ' {output}'
        '&& {conda_wrapper} isoseq3 bcstats'
        ' --json {outdir}/{sample}/{sample}.bcstats.json'
        ' -o {outdir}/{sample}/{sample}.bcstats.tsv'
        ' {output} '
        '&& python3 {basedir}/plot_knees.py'
        ' -t {outdir}/{sample}/{sample}.bcstats.tsv'
        ' -o {outdir}/{sample}/{sample}.cellbackground'
        ' --estimate_percentile 95'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule sort_bam:
    input:
        rules.correct.output
    output:
        "{outdir}/{sample}/{sample}.correct.sort.bam",
    log:
        out='{outdir}/{sample}/5_sort_log.out',
        err='{outdir}/{sample}/5_sort_log.err',
    threads: 20
    resources:
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{conda_wrapper} samtools sort -@ {threads} -t CB'
        ' {input}'
        ' -o {output}'

rule groupdedup:
    input:
        rules.sort_bam.output
    output:
        "{outdir}/{sample}/{sample}.dedup.bam",
    log:
        out='{outdir}/{sample}/6_dedup_log.out',
        err='{outdir}/{sample}/6_dedup_log.err',
    threads: 20
    resources:
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{conda_wrapper} isoseq3 groupdedup'
        ' -j {threads}'
        ' {input}'
        ' {output}'

rule mapping:
    input:
        rules.groupdedup.output
    output:
        "{outdir}/{sample}/{sample}.align.bam",
    log:
        out='{outdir}/{sample}/7_mapping_log.out',
        err='{outdir}/{sample}/7_mapping_log.err',
    threads: 40
    resources:
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{conda_wrapper} pbmm2 align --preset ISOSEQ --sort'
        ' -j {threads}'
        ' {input}'
        ' {genome}'
        ' {output}'

rule collapse:
    input:
        rules.mapping.output
    output:
        "{outdir}/{sample}/{sample}.collapse.gff",
    log:
        out='{outdir}/{sample}/8_collapse_log.out',
        err='{outdir}/{sample}/8_collapse_log.err',
    threads: 20
    resources:
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{conda_wrapper} isoseq3 collapse'
        ' -j {threads}'
        ' {input}'
        ' {output}'

rule sort_gff:
    input:
        rules.collapse.output
    output:
        "{outdir}/{sample}/{sample}.sort.gff",
    log:
        out='{outdir}/{sample}/9_sort_gff_log.out',
        err='{outdir}/{sample}/9_sort_gff_log.err',
    resources:
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{conda_wrapper} pigeon sort'
        ' {outdir}/{sample}/{sample}.collapse.gff'
        ' -o {output}'

rule classify:
    input:
        rules.sort_gff.output
    output:
        "{outdir}/{sample}/{sample}_classification.txt",
    log:
        out='{outdir}/{sample}/10_classify_log.out',
        err='{outdir}/{sample}/10_classify_log.err',
    threads: 20
    params:
        cage=config['cage'],
        polya=config['polya'],
    resources:
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{conda_wrapper} pigeon classify'
        ' -j {threads}'
        ' {input}'
        ' {gtf} {genome}'
        ' --fl {outdir}/{sample}/{sample}.collapse.abundance.txt'
        ' --cage-peak {params.cage}'
        ' --poly-a {params.polya}'
        ' -d {outdir}/{sample} -o {sample}'

rule filter_isoforms:
    input:
        rules.classify.output
    output:
        "{outdir}/{sample}/{sample}_classification.filtered_lite_classification.txt",
    log:
        out='{outdir}/{sample}/11_filter_isoforms_log.out',
        err='{outdir}/{sample}/11_filter_isoforms_log.err',
    threads: 20
    resources:
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{conda_wrapper} pigeon filter'
        ' -j {threads}'
        ' {input}'
        ' --isoforms {outdir}/{sample}/{sample}.sort.gff'

rule report_saturation:
    input:
        rules.filter_isoforms.output
    output:
        "{outdir}/{sample}/{sample}.saturation.txt",
    log:
        out='{outdir}/{sample}/12_report_saturation_log.out',
        err='{outdir}/{sample}/12_report_saturation_log.err',
    threads: 20
    resources:
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{conda_wrapper} pigeon report'
        ' -j {threads}'
        ' {input}'
        ' {output}'


rule make_seurat:
    input:
        rules.filter_isoforms.output
    output:
        "{outdir}/{sample}/genes_seurat/matrix.mtx",
    log:
        out='{outdir}/{sample}/13_make_seurat_log.out',
        err='{outdir}/{sample}/13_make_seurat_log.err',
    threads: 20
    resources:
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{conda_wrapper} pigeon make-seurat'
        ' --keep-ribo-mito-genes --keep-novel-genes'
        ' -j {threads}'
        ' --dedup {outdir}/{sample}/{sample}.dedup.fasta'
        ' --group {outdir}/{sample}/{sample}.collapse.group.txt'
        ' -d {outdir}/{sample}'
        ' -o {sample}'
        ' {input}'

