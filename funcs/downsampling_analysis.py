import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
import argparse
import smm_utils as smm
import plot as pl

parser = argparse.ArgumentParser(description='Runs BMF at downsampled sample-numbers.')
parser.add_argument('-i', '--input',
                    required=True,
                    help='<Required> Input bams to create barcode-umi hashtable.')
parser.add_argument('-o', '--outdir',
                    required=True,
                    help='<Required> Output directory for downsampling files.')
parser.add_argument('-p', '--plotdir',
                    required=True,
                    help='<Required> Output directory to save plots.',
                    )
args = parser.parse_args()

# Seed
np.random.seed(0)

# Make output directory
os.makedirs(args.outdir, exist_ok=True)

# Read input
X = pd.read_csv(args.input, sep='\t', index_col=0)
sample_n = list(range(25,225,25)) + [X.shape[0]]
k_n=[2,3,4,5,6,7,8,9,10]

# Run Downsampling
s_score, rss, evar, kl = smm.downsample_analysis(X, sample_n, n_iter=100, k_n=k_n)

np.save(os.path.join(args.outdir, "s_score.npy"), s_score)
np.save(os.path.join(args.outdir, "rss.npy"), rss)
np.save(os.path.join(args.outdir, "evar.npy"), evar)
np.save(os.path.join(args.outdir, "kl.npy"), kl)

# Heatmap Plots
fig,axes = plt.subplots(2,2,figsize=(16,12))

pl.plot_metric_hm(s_score, sample_n, k_n, ax=axes[0,0], title='Silhouette Score')
pl.plot_metric_hm(rss, sample_n, k_n, ax=axes[0,1], title='Residuals')
pl.plot_metric_hm(evar, sample_n, k_n, ax=axes[1,0], title='Explained Variance')
pl.plot_metric_hm(kl, sample_n, k_n, ax=axes[1,1], title='K-L Divergence')

plt.tight_layout()
plt.savefig(os.path.join(args.plotdir, "figS2a_gridplots_downsample_mean_nmf.pdf"), dpi=100, bbox_inches='tight')

# Metric Factor Grids
pl.plot_dist_per_metric_grid(s_score, sample_n, k_n, y_label='Silhouette Score')
plt.savefig(os.path.join(args.plotdir, "figS2b_sscore_nmf_dist.pdf"), dpi=100, bbox_inches='tight')

pl.plot_dist_per_metric_grid(rss, sample_n, k_n, y_label='Residuals')
plt.savefig(os.path.join(args.plotdir, "figS2b_rss_nmf_dist.pdf"), dpi=100, bbox_inches='tight')

pl.plot_dist_per_metric_grid(evar, sample_n, k_n, y_label='Explained Variance')
plt.savefig(os.path.join(args.plotdir, "figS2b_evar_nmf_dist.pdf"), dpi=100, bbox_inches='tight')

pl.plot_dist_per_metric_grid(kl, sample_n, k_n, y_label='K-L Divergence')
plt.savefig(os.path.join(args.plotdir, "figS2b_skl_nmf_dist.pdf"), dpi=100, bbox_inches='tight')
