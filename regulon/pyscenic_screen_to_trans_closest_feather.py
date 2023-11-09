import pandas as pd
import numpy as np

old_db = "~/genomes/hg38/pyscenic/cisTarget/screen/mc_v10_clust/region_based"
feather =  pd.read_feather(f"{old_db}/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather")


new_db = "/public/home/msu/projects/seq3/20230530/8_regulons/trans/databases"
ids = pd.read_table(f"{new_db}/RIF.hg38_screen.1000bp_up_100bp_down.closest.bed",
                   sep="\t", header=None)
print(ids.iloc[1:10, 1:10])
# extract regions nearest RIF trans and rename region to trans
new_feather = feather.iloc[:, ids[9].values]
new_feather.columns = ids[3].values
print(new_feather.iloc[1:10, 1:10])
# rerank regions
new_feather = new_feather.apply(lambda x: x.rank(method="first").astype('int32'), axis=1)
# add motifs
new_feather["motifs"] = feather["motifs"]
print(new_feather.iloc[1:10, list(range(9)) + [-1]])
new_feather.to_feather(f"{new_db}/RIF_trans_closest.1000bp_up_100bp_down.genes_vs_motifs.rankings.feather")

