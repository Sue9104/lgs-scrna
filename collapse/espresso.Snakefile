import snakemake.utils

snakemake.utils.min_version('6.5.3')

onsuccess:
    print('workflow success')

onerror:
    print('workflow error')

DEFAULT_MEM_MB=4 * 1024  # 4 GB
DEFAULT_TIME_HOURS=12

# Specifying this as an input to a rule will disable that rule.
# This can be used in combination with "ruleorder:" to determine what
# rule should be used to create a particular output file.
UNSATISFIABLE_INPUT='unsatisfiable_input_file_path'

outdir=os.path.abspath(config['outdir'])



def all_input(wildcards):
    inputs = dict()
    inputs.update(run_espresso_q_output())
    if config['enable_visualization']:
        vis_files, vis_dirs = run_visualize_output_files_and_dirs()
        inputs.update(vis_files)
        inputs.update(vis_dirs)

    return inputs


localrules: all
rule all:
    input:
        unpack(all_input),


def get_sample_config_from_wildcards(wildcards):
    sample_name = wildcards.sample
    sample_inputs = config['samples'].get(sample_name)
    return {'name': sample_name, 'inputs': sample_inputs}


def get_sample_input_config_from_wildcards(wildcards):
    sample_config = get_sample_config_from_wildcards(wildcards)
    input_index = int(wildcards.sample_input_i)
    if len(sample_config['inputs']) <= input_index:
        input_config = None
    else:
        input_config = sample_config['inputs'][input_index]

    sample_config['input_i'] = input_index
    sample_config['input'] = input_config
    return sample_config


def get_fast5_dir_path_from_wildcards(wildcards):
    sample_config = get_sample_input_config_from_wildcards(wildcards)
    configured_fast5_dir_path = sample_config['input'].get('fast5_dir')
    if configured_fast5_dir_path is not None:
        return configured_fast5_dir_path

    return UNSATISFIABLE_INPUT


def get_guppy_config_from_wildcards(wildcards):
    sample_config = get_sample_input_config_from_wildcards(wildcards)
    return sample_config['input'].get('guppy_config')


def get_combined_fastq_path_from_wildcards(wildcards):
    sample_config = get_sample_input_config_from_wildcards(wildcards)
    configured_fastq = sample_config['input'].get('fastq')
    if configured_fastq is not None:
        return configured_fastq

    return os.path.join(outdir, 'fastq_dir', sample_config['name'],
                        str(sample_config['input_i']), 'combined.fastq')


def get_sam_path_from_wildcards(wildcards):
    sample_config = get_sample_input_config_from_wildcards(wildcards)
    configured_sam = sample_config['input'].get('sam')
    if configured_sam is not None:
        return configured_sam

    return os.path.join(outdir, 'alignment', sample_config['name'],
                        str(sample_config['input_i']), 'aligned.sam')


def get_bam_path_from_wildcards(wildcards):
    sample_config = get_sample_input_config_from_wildcards(wildcards)
    configured_bam = sample_config['input'].get('bam')
    if configured_bam is not None:
        return configured_bam

    return UNSATISFIABLE_INPUT


def espresso_s_input_sams():
    sams = list()
    # config['samples'] is an unordered dictionary. Sort to be consistent.
    for sample_name in sorted(config['samples']):
        sample_inputs = config['samples'][sample_name]
        for i in range(len(sample_inputs)):
            sam = os.path.join(outdir, 'alignment', sample_name, str(i),
                               'aligned.sorted.sam'),
            sams.append(sam)

    return sams


def reference_file_wildcard_constraints():
    reference_files = config.get('reference_files')
    if reference_files:
        file_names = '|'.join([re.escape(file_name)
                               for file_name in reference_files])
        without_gz = '|'.join([re.escape(file_name[:-3])
                               for file_name in reference_files
                               if file_name.endswith('.gz')])
    else:
        no_match = '^$'  # only matches empty string
        file_names = no_match
        without_gz = no_match

    return {'file_names': file_names, 'without_gz': without_gz}


def get_url_for_download_reference_file(wildcards):
    file_name = wildcards.file_name
    return config['reference_files'][file_name]['url']


rule download_reference_file:
    output:
        ref_file=os.path.join(outdir, 'references', '{file_name}'),
    log:
        out=os.path.join(outdir, 'references',
                         'download_reference_file_{file_name}_log.out'),
        err=os.path.join(outdir, 'references',
                         'download_reference_file_{file_name}_log.err'),
    wildcard_constraints:
        file_name=reference_file_wildcard_constraints()['file_names']
    params:
        url=get_url_for_download_reference_file,
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'curl -L \'{params.url}\''
        ' -o {output.ref_file}'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule unzip_reference_file:
    input:
        gz=os.path.join(outdir, 'references', '{file_name}.gz'),
    output:
        un_gz=os.path.join(outdir, 'references', '{file_name}'),
    log:
        out=os.path.join(outdir, 'references',
                         'unzip_reference_file_{file_name}_log.out'),
        err=os.path.join(outdir, 'references',
                         'unzip_reference_file_{file_name}_log.err'),
    wildcard_constraints:
        file_name=reference_file_wildcard_constraints()['without_gz']
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        ' gunzip -c {input.gz}'
        ' 1> {output.un_gz}'
        ' 2> {log.err}'


def get_base_call_with_guppy_cpu_input_dir_from_wildcards(wildcards):
    if config.get('guppy_gpu'):
        return UNSATISFIABLE_INPUT

    return get_fast5_dir_path_from_wildcards(wildcards)


def get_base_call_with_guppy_gpu_input_dir_from_wildcards(wildcards):
    if not config.get('guppy_gpu'):
        return UNSATISFIABLE_INPUT

    return get_fast5_dir_path_from_wildcards(wildcards)


rule base_call_with_guppy_cpu:
    input:
        fast5_dir=get_base_call_with_guppy_cpu_input_dir_from_wildcards,
    output:
        fastq_dir=directory(os.path.join(outdir, 'fastq_dir', '{sample}',
                                         '{sample_input_i}', 'separate')),
    log:
        out=os.path.join(outdir, 'fastq_dir', '{sample}', '{sample_input_i}',
                         'base_call_with_guppy_cpu_log.out'),
        err=os.path.join(outdir, 'fastq_dir', '{sample}', '{sample_input_i}',
                         'base_call_with_guppy_cpu_log.err'),
    params:
        guppy_config=get_guppy_config_from_wildcards,
        guppy_basecaller_path=os.path.join(config['guppy_bin_path'],
                                           'guppy_basecaller'),
    threads: config['guppy_threads']
    resources:
        mem_mb=config['guppy_mem_gb'] * 1024,
        time_hours=config['guppy_time_hr'],
    shell:
        '{params.guppy_basecaller_path} --input_path {input.fast5_dir}'
        ' --save_path {output.fastq_dir}'
        ' --cpu_threads_per_caller {threads}'
        ' --config {params.guppy_config}'
        ' 1> {log.out}'
        ' 2> {log.err}'


rule base_call_with_guppy_gpu:
    input:
        fast5_dir=get_base_call_with_guppy_gpu_input_dir_from_wildcards,
    output:
        fastq_dir=directory(os.path.join(outdir, 'fastq_dir', '{sample}',
                                         '{sample_input_i}', 'separate')),
    log:
        out=os.path.join(outdir, 'fastq_dir', '{sample}', '{sample_input_i}',
                         'base_call_with_guppy_gpu_log.out'),
        err=os.path.join(outdir, 'fastq_dir', '{sample}', '{sample_input_i}',
                         'base_call_with_guppy_gpu_log.err'),
    params:
        guppy_config=get_guppy_config_from_wildcards,
        guppy_basecaller_path=os.path.join(config['guppy_bin_path'],
                                           'guppy_basecaller'),
    resources:
        mem_mb=config['guppy_mem_gb'] * 1024,
        time_hours=config['guppy_time_hr'],
        gpus=config['guppy_gpus'],
        gpu_name=config['guppy_gpu_name'],
    shell:
        '{params.guppy_basecaller_path} --input_path {input.fast5_dir}'
        ' --save_path {output.fastq_dir}'
        ' --config {params.guppy_config}'
        ' --device au/public/home/msu/projects/seq3/collapse/espresso/results/chrY/q_work_dir/samples.tsv.updated;to'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule combine_fastq_files:
    input:
        fastq_dir=os.path.join(outdir, 'fastq_dir', '{sample}', '{sample_input_i}',
                               'separate'),
    output:
        fastq=os.path.join(outdir, 'fastq_dir', '{sample}', '{sample_input_i}',
                           'combined.fastq'),
    log:
        out=os.path.join(outdir, 'fastq_dir', '{sample}', '{sample_input_i}',
                         'combine_fastq_files_log.out'),
        err=os.path.join(outdir, 'fastq_dir', '{sample}', '{sample_input_i}',
                         'combine_fastq_files_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts', 'combine_fastq_files.py'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} python {params.script}'
        ' --input-dir {input.fastq_dir}'
        ' --out-path {output.fastq}'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule create_bed_from_gtf:
    input:
        gtf=os.path.join(outdir, 'references', config['gtf_name']),
    output:
        bed=os.path.join(outdir, 'references', '{}.bed'.format(config['gtf_name'])),
    log:
        out=os.path.join(outdir, 'references', 'create_bed_from_gtf_log.out'),
        err=os.path.join(outdir, 'references', 'create_bed_from_gtf_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} paftools.js'
        ' gff2bed {input.gtf}'
        ' 1> {output.bed}'
        ' 2> {log.err}'


def align_with_minimap2_input(wildcards):
    inputs = dict()
    inputs['fasta'] = os.path.join(outdir, 'references', config['fasta_name'])
    inputs['fastq'] = get_combined_fastq_path_from_wildcards(wildcards)
    if config['use_annotated_junctions_with_minimap2']:
        inputs['bed'] = os.path.join(
            outdir, 'references', '{}.bed'.format(config['gtf_name']))

    return inputs


def maybe_junc_bed_param(wildcards, input):
    if not config['use_annotated_junctions_with_minimap2']:
        return ''

    return '--junc-bed {}'.format(input.bed)


# If a bam file is configured then use it
ruleorder: convert_bam_to_sam_with_samtools > align_with_minimap2

rule align_with_minimap2:
    input:
        unpack(align_with_minimap2_input),
    output:
        sam=os.path.join(outdir, 'alignment', '{sample}', '{sample_input_i}',
                         'aligned.sam'),
    log:
        out=os.path.join(outdir, 'alignment', '{sample}', '{sample_input_i}',
                         'align_with_minimap2_log.out'),
        err=os.path.join(outdir, 'alignment', '{sample}', '{sample_input_i}',
                         'align_with_minimap2_log.err'),
    benchmark:
        os.path.join(outdir, 'alignment', '{sample}', '{sample_input_i}',
                         'align_with_minimap2_benchmark.out'),
    params:
        conda_wrapper=config['conda_wrapper'],
        preset='splice',
        canonical_splice_sites='b',  # check both strands
        kmer_size=14,
        window_size=4,
        output_secondary_alignments='no',
        maybe_junc_bed=maybe_junc_bed_param,
    threads: config['minimap2_threads']
    resources:
        mem_mb=config['minimap2_mem_gb'] * 1024,
        time_hours=config['minimap2_time_hr'],
    shell:
        '{params.conda_wrapper} minimap2 -a'
        ' -x {params.preset}'
        ' -u{params.canonical_splice_sites}'
        ' -k {params.kmer_size}'
        ' -w {params.window_size}'
        ' {params.maybe_junc_bed}'
        ' --secondary={params.output_secondary_alignments}'
        ' -t {threads}'
        ' {input.fasta} {input.fastq}'
        ' 1> {output.sam}'
        ' 2> {log.err}'

rule convert_bam_to_sam_with_samtools:
    input:
        bam=get_bam_path_from_wildcards,
    output:
        sam=os.path.join(outdir, 'alignment', '{sample}', '{sample_input_i}',
                         'aligned.sam'),
    log:
        out=os.path.join(outdir, 'alignment', '{sample}', '{sample_input_i}',
                         'convert_bam_to_sam_with_samtools_log.out'),
        err=os.path.join(outdir, 'alignment', '{sample}', '{sample_input_i}',
                         'convert_bam_to_sam_with_samtools_log.err'),
    benchmark:
        os.path.join(outdir, 'alignment', '{sample}', '{sample_input_i}',
                         'convert_bam_to_sam_with_samtools_benchmark.out'),
    params:
        conda_wrapper=config['conda_wrapper'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} samtools view'
        ' -h'
        ' -o {output.sam}'
        ' {input.bam}'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule sort_with_samtools:
    input:
        unsorted_sam=get_sam_path_from_wildcards,
    output:
        sorted_sam=os.path.join(outdir, 'alignment', '{sample}', '{sample_input_i}',
                                'aligned.sorted.sam'),
    log:
        out=os.path.join(outdir, 'alignment', '{sample}', '{sample_input_i}',
                         'sort_with_samtools_log.out'),
        err=os.path.join(outdir, 'alignment', '{sample}', '{sample_input_i}',
                         'sort_with_samtools_log.err'),
    benchmark:
        os.path.join(outdir, 'alignment', '{sample}', '{sample_input_i}',
                         'sort_with_samtools_benchmark.out'),
    params:
        conda_wrapper=config['conda_wrapper'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} samtools sort'
        ' -o {output.sorted_sam}'
        ' {input.unsorted_sam}'
        ' 1> {log.out}'
        ' 2> {log.err}'

localrules: write_espresso_s_input_file
rule write_espresso_s_input_file:
    output:
        s_input=os.path.join(outdir, 'run_espresso_s_input.txt'),
    run:
        import os.path
        sams = list()
        names = list()
        # config['samples'] is an unordered dictionary. Sort to be consistent.
        for sample_name in sorted(config['samples']):
            sample_inputs = config['samples'][sample_name]
            names.extend([sample_name] * len(sample_inputs))
            for i in range(len(sample_inputs)):
                sam = os.path.join(outdir, 'alignment', sample_name, str(i),
                                   'aligned.sorted.sam')
                sams.append(sam)

        with open(output[0], 'wt') as handle:
            handle.write('{}\n'.format(','.join(sams)))
            handle.write('{}\n'.format(','.join(names)))

rule run_espresso_s:
    input:
        s_input=os.path.join(outdir, 'run_espresso_s_input.txt'),
        fasta=os.path.join(outdir, 'references', config['fasta_name']),
        gtf=os.path.join(outdir, 'references', config['gtf_name']),
        sams=espresso_s_input_sams(),
    output:
        samples_tsv=os.path.join(outdir, 'samples.tsv'),
        updated_tsv=os.path.join(outdir, 's_work_dir',
                                 'samples.tsv.updated'),
    log:
        out=os.path.join(outdir, 'run_espresso_s_log.out'),
        err=os.path.join(outdir, 'run_espresso_s_log.err'),
    benchmark:
        os.path.join(outdir, 'run_espresso_s_benchmark.out'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join(config['scripts_path'], 'espresso_s_wrapper.py'),
        espresso_s_path=os.path.join(config['espresso_path'], 'ESPRESSO_S.pl'),
        out_dir=lambda wildcards, output: os.path.dirname(output.updated_tsv),
        tmp_dir=outdir,
    threads: config['espresso_s_threads']
    resources:
        mem_mb=config['espresso_s_mem_gb'] * 1024,
        time_hours=config['espresso_s_time_hr'],
    shell:
        '{params.conda_wrapper} python {params.script}'
        ' --s-input {input.s_input}'
        ' --out-tsv {output.samples_tsv}'
        ' --tmp-dir {params.tmp_dir}'
        ' --command'
        ' perl {params.espresso_s_path}'
        ' -A {input.gtf}'
        ' -F {input.fasta}'
        ' -O {params.out_dir}'
        ' -T {threads}'
        ' 1> {log.out}'
        ' 2> {log.err}'

checkpoint split_s_for_c:
    input:
        fasta=os.path.join(outdir, 'references', config['fasta_name']),
        samples_tsv=os.path.join(outdir, 's_work_dir',
                                 'samples.tsv.updated'),
    output:
        samples_tsv=os.path.join(outdir, 'c_work_dir',
                                 'samples.tsv.updated'),
    log:
        out=os.path.join(outdir, 'split_s_for_c_log.out'),
        err=os.path.join(outdir, 'split_s_for_c_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts', 'split_espresso_s_output_for_c.py'),
        orig_work_dir=os.path.join(outdir, 's_work_dir'),
        new_base_dir=os.path.join( outdir, 'c_work_dir'),
        target_reads_per_c=config['target_reads_per_espresso_c_job'],
        sort_memory_buffer_size='2G',
        num_threads_per_c=config['espresso_c_threads'],
    resources:
        mem_mb=config['split_s_for_c_mem_gb'] * 1024,
        time_hours=config['split_s_for_c_time_hr'],
    shell:
        '{params.conda_wrapper} python {params.script}'
        ' --orig-work-dir {params.orig_work_dir}'
        ' --new-base-dir {params.new_base_dir}'
        ' --target-reads-per-c {params.target_reads_per_c}'
        ' --num-threads-per-c {params.num_threads_per_c}'
        ' --genome-fasta {input.fasta}'
        ' --sort-memory-buffer-size {params.sort_memory_buffer_size}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def run_espresso_c_out_dir_param(wildcards, input):
    work_dir = os.path.dirname(input.updated_tsv)
    return os.path.join(work_dir, wildcards.updated_tsv_i)


rule run_espresso_c:
    input:
        updated_tsv=os.path.join(outdir, 'c_work_dir',
                                 'samples.tsv.updated'),
    output:
        done=touch(os.path.join(outdir, 'c_work_dir',
                                'run_espresso_c_{updated_tsv_i}.done')),
    log:
        out=os.path.join(outdir, 'c_work_dir',
                         'run_espresso_c_{updated_tsv_i}_log.out'),
        err=os.path.join(outdir, 'c_work_dir',
                         'run_espresso_c_{updated_tsv_i}_log.err'),
    benchmark:
        os.path.join(outdir, 'c_work_dir',
                     'run_espresso_c_{updated_tsv_i}_benchmark.out'),
    params:
        conda_wrapper=config['conda_wrapper'],
        espresso_c_path=os.path.join(config['espresso_path'], 'ESPRESSO_C.pl'),
        keep_temp='-K' if config['keep_espresso_c_temp'] else '',
        out_dir=run_espresso_c_out_dir_param,
        sub_dir_i='0',
        fasta=os.path.join(outdir, 'c_work_dir', 'fastas',
                           '{updated_tsv_i}.fa'),
    threads: config['espresso_c_threads']
    resources:
        mem_mb=config['espresso_c_mem_gb'] * 1024,
        time_hours=config['espresso_c_time_hr'],
    shell:
        '{params.conda_wrapper}'
        ' perl {params.espresso_c_path}'
        ' -I {params.out_dir}'
        ' -F {params.fasta}'
        ' -X {params.sub_dir_i}'
        ' -T {threads}'
        ' {params.keep_temp}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def combine_c_for_q_input(wildcards):
    inputs = dict()
    split_checkpoint = checkpoints.split_s_for_c.get()
    samples_tsv = split_checkpoint.output.samples_tsv
    with open(samples_tsv, 'rt') as handle:
        for line in handle:
            columns = line.strip().split('\t')
            c_index = columns[2]
            inputs[c_index] = os.path.join(
                outdir, 'c_work_dir',
                'run_espresso_c_{}.done'.format(c_index))

    return inputs


rule combine_c_for_q:
    input:
        unpack(combine_c_for_q_input),
    output:
        updated_tsv=os.path.join(outdir, 'q_work_dir',
                                 'samples.tsv.updated'),
    log:
        out=os.path.join(outdir, 'combine_c_for_q_log.out'),
        err=os.path.join(outdir, 'combine_c_for_q_log.err'),
    benchmark:
        os.path.join(outdir, 'combine_c_for_q_benchmark.out'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts', 'combine_espresso_c_output_for_q.py'),
        c_work_dir=os.path.join(outdir, 'c_work_dir'),
        q_work_dir=os.path.join(outdir, 'q_work_dir'),
    resources:
        mem_mb=config['combine_c_for_q_mem_gb'] * 1024,
        time_hours=config['combine_c_for_q_time_hr'],
    shell:
        '{params.conda_wrapper} python {params.script}'
        ' --c-base-work-dir {params.c_work_dir}'
        ' --new-base-dir {params.q_work_dir}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def run_espresso_q_output():
    outputs = dict()
    outputs['abundance_esp'] = os.path.join(outdir, 'q_work_dir',
                                            'samples_N2_R0_abundance.esp')
    outputs['updated_gtf'] = os.path.join(outdir, 'q_work_dir',
                                          'samples_N2_R0_updated.gtf')
    if config['output_compatible_isoforms']:
        outputs['isoform_tsv'] = os.path.join(
            outdir, 'q_work_dir', 'samples_N2_R0_compatible_isoform.tsv')

    return outputs


def maybe_compat_iso_param(wildcards, output):
    if not config['output_compatible_isoforms']:
        return ''

    return '-V {}'.format(output.isoform_tsv)


rule run_espresso_q:
    input:
        updated_tsv=os.path.join(outdir, 'q_work_dir',
                                 'samples.tsv.updated'),
        gtf=os.path.join(outdir, 'references', config['gtf_name']),
    output:
        **run_espresso_q_output()
    log:
        out=os.path.join(outdir, 'run_espresso_q_log.out'),
        err=os.path.join(outdir, 'run_espresso_q_log.err'),
    benchmark:
        os.path.join(outdir, 'run_espresso_q_benchmark.out'),
    params:
        conda_wrapper=config['conda_wrapper'],
        espresso_q_path=os.path.join(config['espresso_path'], 'ESPRESSO_Q.pl'),
        maybe_compat_iso=maybe_compat_iso_param,
    threads: config['espresso_q_threads']
    resources:
        mem_mb=config['espresso_q_mem_gb'] * 1024,
        time_hours=config['espresso_q_time_hr'],
    shell:
        '{params.conda_wrapper} perl {params.espresso_q_path}'
        ' -A {input.gtf}'
        ' -L {input.updated_tsv}'
        ' -T {threads}'
        ' {params.maybe_compat_iso}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def run_visualize_output_files_and_dirs():
    files = dict()
    dirs = dict()

    files['fasta_index'] = os.path.join(outdir, 'references',
                                        '{}.fai'.format(config['fasta_name']))
    base_path = 'visualization'
    dirs['target_genes'] = os.path.join(base_path, 'target_genes')
    for sample_name in config['samples']:
        files['{}_bigwig'.format(sample_name)] = os.path.join(
            base_path, '{}.bw'.format(sample_name))

    return files, dirs


def run_visualize_output():
    files, dirs = run_visualize_output_files_and_dirs()
    for d_name, d_path in dirs.items():
        files[d_name] = directory(d_path)

    return files


def visualize_normalize_param():
    if config['vis_normalize_counts']:
        return '--normalize-counts-to-cpm'

    return ''


def visualize_output_dir(wildcards, output):
    first_sample_name = None
    for sample_name in config['samples']:
        first_sample_name = sample_name

    output_file_name = '{}_bigwig'.format(first_sample_name)
    output_file_path = output.get(output_file_name)
    return os.path.dirname(output_file_path)


rule run_visualize:
    input:
        fasta=os.path.join(outdir, 'references', config['fasta_name']),
        updated_gtf=os.path.join(outdir, 'q_work_dir',
                                 'samples_N2_R0_updated.gtf'),
        abundance_esp=os.path.join(outdir, 'q_work_dir',
                                   'samples_N2_R0_abundance.esp'),
    output:
        **run_visualize_output()
    log:
        out=os.path.join(outdir, 'visualization', 'run_visualize_log.out'),
        err=os.path.join(outdir, 'visualization', 'run_visualize_log.err'),
    benchmark:
        os.path.join(outdir, 'visualization', 'run_visualize_benchmark.out'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join(config['visualization_path'], 'visualize.py'),
        target_gene=config['vis_target_gene'],
        minimum_count=config['vis_minimum_count'],
        descriptive_name=config['vis_descriptive_name'],
        normalize=visualize_normalize_param(),
        output_dir=visualize_output_dir,
    resources:
        mem_mb=config['visualize_mem_gb'] * 1024,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} python {params.script}'
        ' --genome-fasta {input.fasta}'
        ' --updated-gtf {input.updated_gtf}'
        ' --abundance-esp {input.abundance_esp}'
        ' --target-gene {params.target_gene}'
        ' --minimum-count {params.minimum_count}'
        ' --descriptive-name {params.descriptive_name}'
        ' --output-dir {params.output_dir}'
        ' {params.normalize}'
        ' 1> {log.out}'
        ' 2> {log.err}'
