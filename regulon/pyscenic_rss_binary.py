import argparse
import shutil
parser = argparse.ArgumentParser( description='Read pyscenic indir and export RSS and Binary results')
parser.add_argument( 'indir', type=str, help='input  dir')
parser.add_argument( 'outdir', type=str, help='output  dir')
parser.add_argument( 'regulons', type=str, help='regulons csv')
parser.add_argument( 'aucell', type=str, help='aucell csv')
args = parser.parse_args()
print(args.indir)


import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
import operator as op
from cytoolz import compose

import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss


indir = args.indir
outdir = args.outdir
os.makedirs(outdir, exist_ok=True)
os.chdir( outdir )


# enriched regulons
def derive_regulons(motifs):
    #motifs.columns = motifs.columns.droplevel(0)

    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f

    # For the creation of regulons we only keep the 10-species databases and the activating modules. We also remove the
    # enriched motifs for the modules that were created using the method 'weight>50.0%' (because these modules are not part
    # of the default settings of modules_from_adjacencies anymore.
    print(motifs['Context'].head())
    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains('weight>50.0%')), motifs.Context), dtype=bool) & \
        #np.fromiter(map(contains(*db_names), motifs.Context), dtype=np.bool) & \
        np.fromiter(map(contains('activating'), motifs['Context']), dtype=bool)]

    # We build regulons only using enriched motifs with a NES of 3.0 or higher; we take only directly annotated TFs or TF annotated
    # for an orthologous gene into account; and we only keep regulons with at least 10 genes.
    regulons = list(filter(
        lambda r: len(r) >= 10,
            df2regulons(motifs[(motifs['NES'] >= 3.0)
            & ((motifs['Annotation'] == 'gene is directly annotated')
            | (motifs['Annotation'].str.startswith('gene is orthologous to')
            & motifs['Annotation'].str.endswith('which is directly annotated for motif')))])))
    # Rename regulons, i.e. remove suffix.
    return list(map(lambda r: r.rename(r.transcription_factor), regulons))

MOTIFS_FNAME = f"{args.indir}/{args.regulons}"
shutil.copy2(MOTIFS_FNAME, "regulons.csv")
df_motifs = load_motifs(MOTIFS_FNAME)
print(df_motifs.head())
#regulons = derive_regulons(df_motifs)
#REGULONS_DAT_FNAME = "regulons.dat"
#with open(REGULONS_DAT_FNAME, 'wb') as f:
#    pickle.dump(regulons, f)


shutil.copy2(f"{indir}/celltype.genes_fine.tsv", "celltype.genes_fine.tsv")
shutil.copy2(f"{indir}/celltype.trans_fine.tsv", "celltype.trans_fine.tsv")
cellAnnot_genes = pd.read_table("celltype.genes_fine.tsv", sep='\t', index_col=0)
cellAnnot_trans = pd.read_table("celltype.trans_fine.tsv", sep='\t', index_col=0)
# need to transpose when using pipe cmd shell
auc_mtx = pd.read_csv(f"{args.indir}/{args.aucell}", index_col=0).transpose()
auc_mtx.to_csv("auc_mtx.normal.csv")
#auc_mtx = pd.read_csv("auc_mtx.short.csv", index_col=0)

# Regulon specificity scores (RSS) across predicted cell types
rss_genes = regulon_specificity_scores(auc_mtx, cellAnnot_genes["cell_type"])
rss_genes.to_csv("genes-cellType-RSS.csv")
print(rss_genes.head())
rss_trans = regulon_specificity_scores(auc_mtx, cellAnnot_trans["cell_type"])
rss_trans.to_csv("trans-cellType-RSS.csv")

# calculate z-score for each regulon
auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
for col in list(auc_mtx.columns):
    auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)
auc_mtx_Z.to_csv("auc_mtx_Z.csv")

# Regulon specificity scores (RSS) across predicted cell types
print("binarizing...")
#bin_mtx, thresholds = binarize(auc_mtx, num_workers=20)
bin_mtx, thresholds = binarize(auc_mtx)
BIN_MTX_FNAME = "bin.csv"
THR_FNAME = "thresholds.csv"
bin_mtx.to_csv(BIN_MTX_FNAME)
thresholds.to_frame().rename(columns={0:'threshold'}).to_csv(THR_FNAME)
## futher plot for specific regulons
#bin_mtx = pd.read_csv(BIN_MTX_FNAME, index_col=0)
#thresholds = pd.read_csv(THR_FNAME, index_col=0).threshold
#fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, figsize=(8, 4), dpi=100)
#
#plot_binarization(auc_mtx, 'NFKB2', thresholds['NFKB2'], ax=ax1)
#plot_binarization(auc_mtx, 'MITF', thresholds['MITF'], ax=ax2)
#plot_binarization(auc_mtx, 'FOXP3', thresholds['FOXP3'], ax=ax3)
#plot_binarization(auc_mtx, 'PAX5', thresholds['PAX5'], ax=ax4)
#plot_binarization(auc_mtx, 'IRF8', thresholds['IRF8'], ax=ax5)
#plot_binarization(auc_mtx, 'IRF3', thresholds['IRF3'], ax=ax6)
#plot_binarization(auc_mtx, 'MLX', thresholds['MLX'], ax=ax7)
#plot_binarization(auc_mtx, 'YY1', thresholds['YY1'], ax=ax8)
#
#plt.tight_layout()
#savesvg('hists - GSE115978 - binarization.svg', fig)



