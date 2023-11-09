database=/public/home/msu/genomes/hg38/pyscenic/
softdir=/public/home/msu/softwares/pyscenic
cores=20

indir=$(dirname $1)
infile=$(basename $1)
echo ${indir}
echo ${infile}

echo singularity run --bind ${database}:/database,${indir}:/data ${softdir}/aertslab-pyscenic-0.12.1.sif pyscenic grn -t -o /data/expr_mat.adjacencies.tsv /data/${infile} /database/TFlist/allTFs_hg38.txt --num_workers ${cores}
singularity run --bind ${database}:/database,${indir}:/data ${softdir}/aertslab-pyscenic-0.12.1.sif pyscenic grn -t -o /data/expr_mat.adjacencies.tsv /data/${infile} /database/TFlist/allTFs_hg38.txt --num_workers ${cores}

echo singularity run --bind ${database}:/database,${indir}:/data ${softdir}/aertslab-pyscenic-0.12.1.sif pyscenic ctx -t /data/expr_mat.adjacencies.tsv /database/cisTarget/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /database/cisTarget/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname /database/motif2TF/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /data/${infile} --output /data/regulons.csv --num_workers ${cores} -a
singularity run --bind ${database}:/database,${indir}:/data ${softdir}/aertslab-pyscenic-0.12.1.sif pyscenic ctx -t /data/expr_mat.adjacencies.tsv /database/cisTarget/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /database/cisTarget/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname /database/motif2TF/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /data/${infile} --output /data/regulons.csv --num_workers ${cores} -a

echo singularity run --bind ${indir}:/data ${softdir}/aertslab-pyscenic-0.12.1.sif pyscenic aucell -t /data/${infile} /data/regulons.csv -o /data/auc_mtx.csv --num_workers ${cores}
singularity run --bind ${indir}:/data ${softdir}/aertslab-pyscenic-0.12.1.sif pyscenic aucell -t /data/${infile} /data/regulons.csv -o /data/auc_mtx.csv --num_workers ${cores}
