basedir=$(dirname "$0")
intsv=$1
outdir=$2
#genome=$3
#gtf=$4
genome=/public/home/msu/genomes/hg38/ensembl/Homo_sapiens.GRCh38.dna_sm.primary_assembly.chr.fa
gtf=/public/home/msu/genomes/hg38/gencode/gencode.v42.primary_assembly.annotation.gtf
bin=/public/home/msu/softwares/espresso
cores=10
memory=32G
# create outdir if not exists
[ ! -d "${outdir}" ] && mkdir ${outdir}

# preprocess
echo perl ${bin}/src/ESPRESSO_S.pl -L ${intsv} -F ${genome} -A ${gtf} -O ${outdir}/s_work_dir
perl ${bin}/src/ESPRESSO_S.pl -L ${intsv} -F ${genome} -A ${gtf} -O ${outdir}/s_work_dir

# split for collapse
echo python ${bin}/snakemake/scripts/split_espresso_s_output_for_c.py --orig-work-dir ${outdir}/s_work_dir --new-base-dir ${outdir}/c_work_dir --target-reads-per-c 1000 --genome-fasta ${genome} --sort-memory-buffer-size ${memory}
\rm -rf ${outdir}/c_work_dir && python ${bin}/snakemake/scripts/split_espresso_s_output_for_c.py --orig-work-dir ${outdir}/s_work_dir --new-base-dir ${outdir}/c_work_dir --target-reads-per-c 1000 --genome-fasta ${genome} --sort-memory-buffer-size ${memory}

# collapse
samples=`wc -l ${intsv}| awk '{print $1}'`
loops=$(($samples<5?$samples:5))
for (( i=0; i < ${loops}; ++i ))
do
  echo perl ${bin}/src/ESPRESSO_C.pl -I ${outdir}/c_work_dir/$i -F ${outdir}/c_work_dir/fastas/$i.fa -X 0 -T ${cores}
  perl ${bin}/src/ESPRESSO_C.pl -I ${outdir}/c_work_dir/$i -F ${outdir}/c_work_dir/fastas/$i.fa -X 0 -T ${cores}
done

# combine for quantify
echo python ${bin}/snakemake/scripts/combine_espresso_c_output_for_q.py --c-base-work-dir ${outdir}/c_work_dir --new-base-dir ${outdir}/q_work_dir
python ${bin}/snakemake/scripts/combine_espresso_c_output_for_q.py --c-base-work-dir ${outdir}/c_work_dir --new-base-dir ${outdir}/q_work_dir

# quantify
echo perl ${bin}/src/ESPRESSO_Q.pl -A ${gtf} -L ${outdir}/q_work_dir/samples.tsv.updated -V ${outdir}/q_work_dir/compatible_isoform.tsv
perl ${bin}/src/ESPRESSO_Q.pl -A ${gtf} -L ${outdir}/q_work_dir/samples.tsv.updated -V ${outdir}/q_work_dir/compatible_isoform.tsv
