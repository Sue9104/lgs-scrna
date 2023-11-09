import pandas as pd
import numpy as np

feather =  pd.read_feather("~/genomes/hg38/pyscenic/cisTarget/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather")
ids = pd.read_table("/public/home/msu/projects/seq3/20230530/8_regulons/trans/databases/RIF.hg38_screen.1000bp_up_100bp_down.overlap_40.bed",
                   sep="\t", header=None)
ids

new_feather = feather.iloc[:, ids[3].values]
new_feather.columns = ids[9].values
new_feather = new_feather.apply(lambda x: x.rank(method="first").astype('int32'), axis=1)
new_feather["motifs"] = feather["motifs"]

new_feather.to_feather("/public/home/msu/genomes/hg38/pyscenic/cisTarget/screen/RIF.trans.hg38_screen_v10_clust.1000bp_up_100bp_down.feather")

