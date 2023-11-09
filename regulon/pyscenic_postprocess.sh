database=/public/home/msu/genomes/hg38/pyscenic/
softdir=/public/home/msu/softwares/pyscenic
pipedir=/public/home/msu/pipelines/lgs-scrna/regulon
cores=20

indir=$(realpath -s $1)
echo ${indir}
outdir=$2
regulons=$(basename $3)
aucell=$(basename $4)

echo singularity run -B ${pipedir}:/pipelines,${indir}:/data /public/home/msu/softwares/pyscenic/aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif python //pipelines/pyscenic_rss_binary.py /data ${outdir} ${regulons} ${aucell}
singularity run -B ${pipedir}:/pipelines,${indir}:/data /public/home/msu/softwares/pyscenic/aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif python //pipelines/pyscenic_rss_binary.py /data ${outdir} ${regulons} ${aucell}
