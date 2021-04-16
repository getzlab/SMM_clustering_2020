import pandas as pd
import sys
import os
import numpy as np
import signatureanalyzer as sa
from typing import Union
import nimfa
from tqdm import tqdm
import sklearn
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import warnings
warnings.filterwarnings("ignore")

def bnmf(df: pd.DataFrame, K: int = 4, n_run: int = 10, **kwargs):
    """
    Binary matrix factorization wrapper.
    ----------------------
    Args:
        * pd.DataFrame: (features x samples)

    Returns:
        * H: pd.Dataframe (n_samples x K)
        * W: pd.DataFrame (K x N_features)
    """
    bmf = nimfa.Bmf(df.values, rank=K, n_run=n_run, **kwargs)
    bmf_fit = bmf()

    W = pd.DataFrame(bmf_fit.fit.W, index=df.index).T
    H = pd.DataFrame(bmf_fit.fit.H, columns=df.columns).T

    H.columns = H.columns.astype(str)
    W.index = W.index.astype(str)

    W,H = sa.utils.select_signatures(W.T,H.T)

    return H, W, bmf, bmf_fit

def consensus_cluster(H_matrices: list):
    """
    Consensus clustering of bnmf results.
    -----------------------
    Args:
        * filepath: path to the output .5 file of ARD-NMF runs
    Returns:
        * pd.DataFrame: consensus matrix from results
        * pd.Series: assignment probability for selected cluster
    """

    x = np.vstack([df.loc[:,'max_id'].values for df in H_matrices])
    consensus_matrix = np.vstack([(x[:,[y]] == x[:]).sum(0) for y in range(x.shape[1])])

    df = pd.DataFrame(consensus_matrix, index=H_matrices[0].index, columns=H_matrices[0].index)
    df = df.loc[H_matrices[0].sort_values('max_id').index, H_matrices[0].sort_values('max_id').index]

    assign_p = pd.concat([df.loc[
        H_matrices[0][H_matrices[0]['max_id']==x].index,
        H_matrices[0][H_matrices[0]['max_id']==x].index
    ].mean(1) for x in set(H_matrices[0]['max_id'])])

    assign_p.name = 'assignment'
    return df, assign_p

def fisher_exact(
    X: pd.DataFrame,
    metadata: pd.DataFrame,
    groupby: str,
    fdr_alpha: float = 0.05,
    fdr_method: str = 'fdr_bh',
    **kwargs
    ):
    """
    Fisher Exact Test.
    -------------------
    Performs fisher exact test by comparing proportions of binary features for full
    population and cluster specific populations:

        (present vs absent)
        (within cluster vs outside cluster)

    Args:
        * X: input binary matrix
        * metadata: metadata for sample set
        * groupby: clustering to compute exact test for
        * fdr_alpha: FDR correction thresh
        * fdr_method: FDR method (statsmodels)
        ** kwargs: for exact test (stats.fisher_exact)

    Returns:
        * pd.DataFrame: results with pval, adj_pval, and oddsratio
    """
    from statsmodels.stats.multitest import multipletests
    import scipy.stats as stats

    # Present & Absent
    X_present = X.sum(0)
    X_absent = X.shape[0] - X.sum(0)

    # Within cluster & Out cluster
    X_ci = X.join(metadata).groupby(groupby).sum()[X.columns]
    X_co = X_present - X_ci

    # Initialize results
    pval = X_ci.T.copy()
    pval_adj = X_ci.T.copy()
    odds_r = X_ci.T.copy()

    # Perform fisher exact test
    for clust in X_ci.index:
        for feat in X_ci.columns:
            odds_r.loc[feat,clust],pval.loc[feat,clust] = stats.fisher_exact([[X_present[feat], X_absent[feat]], [X_ci.loc[clust,feat], X_co.loc[clust,feat]]], alternative='less', **kwargs)

        _,pval_adj[clust],_,_ = multipletests(pval[clust], alpha=fdr_alpha, method=fdr_method)

    # Melt results
    pval_adj = pd.melt(pval_adj.reset_index(), id_vars=['index'], value_vars=pval_adj.columns).rename(
        columns={'index':'feat', 'value':'pval_adj'}).set_index(['feat',groupby])

    pval = pd.melt(pval.reset_index(), id_vars=['index'], value_vars=pval.columns).rename(
        columns={'index':'feat', 'value':'pval'}).set_index(['feat',groupby])

    odds_r = pd.melt(odds_r.reset_index(), id_vars=['index'], value_vars=odds_r.columns).rename(
        columns={'index':'feat', 'value':'odds_r'}).set_index(['feat',groupby])

    return pval.join(pval_adj).join(odds_r)

def downsample_analysis(
    X: pd.DataFrame,
    sample_n: list,
    n_iter: int = 100,
    k_n: list = list(range(2,11)),
    seed=100
    ):
    """
    Downsampling analysis.
    ----------------
    Args:
        * X: pd.DataFrame (samples x features)
        * sample_n: list of integers (n samples to downsample)
        * n_iter: number of bootstrap

    Returns:
        * tuple:
            silhoutte score (dice similarity)
            rss score
            evar
            kl divergence
    """
    from sklearn.metrics import silhouette_score

    if seed is not None:
        np.random.seed(seed)

    samples = np.array(X.index)
    s_score = np.zeros((len(sample_n), len(k_n), n_iter))
    rss_score = np.zeros((len(sample_n), len(k_n), n_iter))
    evar_score = np.zeros((len(sample_n), len(k_n), n_iter))
    kl_score = np.zeros((len(sample_n), len(k_n), n_iter))

    for s_idx,s in enumerate(sample_n):
        for i in tqdm(range(n_iter), desc="n = {}".format(s)):
            idx = np.random.choice(samples, s, replace=False)
            X_run = X.loc[idx].T

            for k_idx,k in enumerate(k_n):
                H, W, bmf, bmf_fit = bnmf(X_run, K=k, n_run=1, seed=None)
                s_score[s_idx,k_idx,i] = sklearn.metrics.silhouette_score(X_run.T.astype('boolean'), H['max_id'], metric='dice')
                rss_score[s_idx,k_idx,i] = bmf_fit.summary()['rss']
                evar_score[s_idx,k_idx,i] = bmf_fit.summary()['evar']
                kl_score[s_idx,k_idx,i] = bmf_fit.summary()['kl']

    return s_score,rss_score,evar_score,kl_score

#------------------------------------------------------------------
# RNA Helpers
#------------------------------------------------------------------
def tpm_loader(tpm, counts, samples=None, filter_thresh=True):
    """
    Bulk load dataset.
    """
    from qtl.norm import deseq2_size_factors

    # Load data
    tpm = pd.read_csv(tpm, sep='\t', skiprows=2, index_col=0)
    counts = pd.read_csv(counts, sep='\t', skiprows=2, index_col=0)
    gene_name = tpm.loc[:,['Description']]
    tpm = tpm.iloc[:,1:]

    if samples is not None:
        tpm = tpm.loc[:,samples]

    # Filter counts
    if filter_thresh:
        tpm = tpm[(np.sum(tpm >= 0.1, 1) > tpm.shape[1]*0.2) & (np.sum(counts.iloc[:,1:] >= 6, 1) > tpm.shape[1]*0.2)]

    return tpm, np.log2(1+tpm / deseq2_size_factors(tpm)), counts, gene_name

#------------------------------------------------------------------
# From Francois
#------------------------------------------------------------------
def get_pcs(gct_df, normalize=True, C=None, n_components=5, return_genes=False):
    """
    Scale input GCT, threshold, normalize and calculate PCs
    """
    if normalize:
        gct_norm_std_df = normalize_counts(gct_df, C=C)
    else:
        gct_norm_std_df = gct_df

    pca = sklearn.decomposition.PCA(n_components=n_components)
    pca.fit(gct_norm_std_df.T)
    P = pca.transform(gct_norm_std_df.T)
    P_df = pd.DataFrame(P, index=gct_norm_std_df.columns)

    if return_genes:
        return P_df, pca, gct_norm_std_df.index.values
    else:
        return P_df, pca

def plot_pca(P_df, pca, c=None, cohort_s=None, cohort_colors=None, cohort_args=None, order=[1,2,3], outliers=None, title='',
    vmin=None, vmax=None, alpha=1, lw=0, s=30, cmap=plt.cm.Spectral_r, cticks=None, cticklabels=None, clabel='',
    show_legend=True, show_ax2=True):
    """
    cohort_s: Series encoding cohorts
    cohort_colors: dict

    Modes:
    """
    if cohort_s is not None:
        cohorts = cohort_s.unique()
        nc = len(cohorts)
        if cohort_colors is None and cohort_args is None:
            # cohort_colors = {i:j for i,j in zip(cohorts, cm.get_cmap(cmap, nc)(np.arange(nc)))}
            cohort_colors = {i:j for i,j in zip(cohorts, sns.husl_palette(nc, s=1, l=0.6))}
        if cohort_args is None:
            cohort_args = {}
            for k in np.unique(cohort_s):
                cohort_args[k] = {'color': cohort_colors[k], 'marker':'o', 'edgecolor':'none', 's':s}

    if show_ax2:
        fig = plt.figure(facecolor=(1,1,1), figsize=(10.5,5.5))
        ax1 = fig.add_axes(np.array([1/10.5, 0.75/5.5, 4/10.5, 4/5.5]))
    else:
        fig = plt.figure(facecolor=(1,1,1), figsize=(5.5,5.5))
        ax1 = fig.add_axes(np.array([1/5.5, 0.75/5.5, 4/5.5, 4/5.5]))
    if cohort_s is None:  # c[P_df.index]
        sa = ax1.scatter(P_df[order[1]-1], P_df[order[0]-1], c=c, cmap=cmap, vmin=vmin, vmax=vmax, lw=lw, alpha=alpha, s=s)
    else:
        for k in np.unique(cohort_s):
        # for k in cohort_s.unique():
            i = cohort_s[cohort_s==k].index
            ax1.scatter(P_df.loc[i,order[1]-1], P_df.loc[i,order[0]-1], alpha=alpha, label=k, **cohort_args[k])
    format_plot(ax1, fontsize=10)
    ax1.set_xlabel('PC {0} ({1:.2f}%)'.format(order[1], pca.explained_variance_ratio_[order[1]-1]*100), fontsize=12)
    ax1.set_ylabel('PC {0} ({1:.2f}%)'.format(order[0], pca.explained_variance_ratio_[order[0]-1]*100), fontsize=12)

    if show_ax2:
        ax2 = fig.add_axes(np.array([6/10.5, 0.75/5.5, 4/10.5, 4/5.5]))
        if cohort_s is None:
            ax2.scatter(P_df[order[2]-1], P_df[order[0]-1], c=c, cmap=cmap, vmin=vmin, vmax=vmax, lw=lw, alpha=alpha, s=s)
        else:
            for k in np.unique(cohort_s):
                i = cohort_s[cohort_s==k].index
                ax2.scatter(P_df.loc[i,order[2]-1], P_df.loc[i,order[0]-1], alpha=alpha, label=k, **cohort_args[k])
            # ax2.legend(loc=3, fontsize=10, scatterpoints=1, handletextpad=0.1, framealpha=0.5, bbox_to_anchor=(-0.5,-0.1))

        format_plot(ax2, fontsize=10)
        ax2.set_xlabel('PC {0} ({1:.2f}%)'.format(order[2], pca.explained_variance_ratio_[order[2]-1]*100), fontsize=12)
        ax2.set_ylabel('PC {0} ({1:.2f}%)'.format(order[0], pca.explained_variance_ratio_[order[0]-1]*100), fontsize=12)

    if outliers is not None:
        ax1.scatter(P_df.loc[outliers, order[1]-1], P_df.loc[outliers, order[0]-1], c='none', edgecolors='r', marker='s', lw=1, alpha=1, s=50, label=None)
        if show_ax2:
            ax2.scatter(P_df.loc[outliers, order[2]-1], P_df.loc[outliers, order[0]-1], c='none', edgecolors='r', marker='s', lw=1, alpha=1, s=50, label=None)

    fig.suptitle(title, fontsize=12)

    if cohort_s is not None and show_legend:
        # ax2.legend(loc=0, fontsize=10, scatterpoints=1, handletextpad=0.1, framealpha=0.5, bbox_to_anchor=(-0.5,-0.1))
        leg = ax1.legend(loc=0, fontsize=9, scatterpoints=1, handletextpad=0.1, framealpha=1, labelspacing=0.35)
        for lh in leg.legendHandles:
            lh.set_alpha(1)

    # if cohort_s is None and c is not None and not isinstance(c, list) and not isinstance(c, str):
    if cohort_s is None and c is not None and len(c)==P_df.shape[0]:
        if show_ax2:
            cax = fig.add_axes(np.array([3.5/10.5, 5/5.5, 1.5/10.5, 0.15/5.5]))
        else:
            cax = fig.add_axes(np.array([3.5/5.5, 5/5.5, 1.5/5.5, 0.15/5.5]))
        # cax = fig.add_axes(np.array([3.5/10.5, 4.85/5.5, 1.5/10.5, 0.15/5.5]))
        hc = plt.colorbar(sa, cax=cax, orientation='horizontal')
        if cticks is not None:
            hc.set_ticks(cticks)
        if cticklabels is not None:
            # hc.set_ticks([0,0.5,1])
            hc.ax.tick_params(labelsize=9)
            # cax.invert_xaxis()
            cax.set_xticklabels(cticklabels, fontsize=10)

        hc.locator = ticker.MaxNLocator(integer=True, min_n_ticks=2, nbins=5)
        hc.update_ticks()

        cax.set_ylabel(clabel, rotation=0, ha='right', va='center', fontsize=12)
    return fig

def format_plot(ax, tick_direction='out', tick_length=4, hide=['top', 'right'], hide_spines=True, lw=1, fontsize=9):

    for i in ['left', 'bottom', 'right', 'top']:
        ax.spines[i].set_linewidth(lw)

    # ax.axis["left"].major_ticklabels.set_ha("left")
    ax.tick_params(axis='both', which='both', direction=tick_direction, labelsize=fontsize)

    # set tick positions
    if 'top' in hide and 'bottom' in hide:
        ax.get_xaxis().set_ticks_position('none')
    elif 'top' in hide:
        ax.get_xaxis().set_ticks_position('bottom')
    elif 'bottom' in hide:
        ax.get_xaxis().set_ticks_position('top')
    else:
        ax.get_xaxis().set_ticks_position('both')

    if 'left' in hide and 'right' in hide:
        ax.get_yaxis().set_ticks_position('none')
    elif 'left' in hide:
        ax.get_yaxis().set_ticks_position('right')
    elif 'right' in hide:
        ax.get_yaxis().set_ticks_position('left')
    else:
        ax.get_yaxis().set_ticks_position('both')

    if hide_spines:
        for i in hide:
            ax.spines[i].set_visible(False)


    # adjust tick size
    for line in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
    #for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(tick_length) # tick length
        line.set_markeredgewidth(lw) # tick line width

    for line in (ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True)):
        line.set_markersize(tick_length/2) # tick length
        line.set_markeredgewidth(lw/2) # tick line width

def plot_pca_ax(
    P_df,
    pca,
    ax=None,
    c=None,
    cohort_s=None,
    cohort_colors=None,
    cohort_args=None,
    order=[1,2,3],
    outliers=None,
    title='',
    vmin=None,
    vmax=None,
    alpha=1,
    lw=0,
    s=30,
    cmap=plt.cm.Spectral_r,
    cticks=None,
    cticklabels=None,
    clabel='',
    show_legend=True,
    plot_color_bar=True
    ):
    """
    PCA Plot by axis.
    -------------------
    cohort_s: Series encoding cohorts
    cohort_colors: dict
    Modes:
    """
    if cohort_s is not None:
        cohorts = cohort_s.unique()
        nc = len(cohorts)
        if cohort_colors is None and cohort_args is None:
            cohort_colors = {i:j for i,j in zip(cohorts, sns.husl_palette(nc, s=1, l=0.6))}
        if cohort_args is None:
            cohort_args = {}
            for k in np.unique(cohort_s):
                cohort_args[k] = {'color': cohort_colors[k], 'marker':'o', 'edgecolor':'none', 's':s}

    if ax is None:
        fig,ax = plt.subplots(figsize=(6,6))

    if cohort_s is None:
        sa = ax.scatter(P_df[order[1]-1], P_df[order[0]-1], c=c, cmap=cmap, vmin=vmin, vmax=vmax, lw=lw, alpha=alpha, s=s)
    else:
        for k in np.unique(cohort_s):
            i = cohort_s[cohort_s==k].index
            ax.scatter(P_df.loc[i,order[1]-1], P_df.loc[i,order[0]-1], alpha=alpha, label=k, **cohort_args[k])

    format_plot(ax, fontsize=10)
    ax.set_xlabel('PC {0} ({1:.2f}%)'.format(order[1], pca.explained_variance_ratio_[order[1]-1]*100), fontsize=12)
    ax.set_ylabel('PC {0} ({1:.2f}%)'.format(order[0], pca.explained_variance_ratio_[order[0]-1]*100), fontsize=12)

    if outliers is not None:
        ax.scatter(P_df.loc[outliers, order[1]-1], P_df.loc[outliers, order[0]-1], c='none', edgecolors='r', marker='s', lw=1, alpha=1, s=50, label=None)

    ax.set_title(title, fontsize=12)

    if cohort_s is not None and show_legend:
        leg = ax.legend(loc=0, fontsize=6, scatterpoints=1, handletextpad=0.1, framealpha=1, labelspacing=0.35)
        for lh in leg.legendHandles:
            lh.set_alpha(1)

    if cohort_s is None and c is not None and len(c)==P_df.shape[0]:
        x1 = ax.get_position().x1
        y1 = ax.get_position().y1

        if plot_color_bar:
            fig = plt.gcf()
            cax = fig.add_axes(np.array([x1, y1*5/5.5, 0.15/5.5, 1/5.5]))
            hc = plt.colorbar(sa, cax=cax, orientation='vertical')

            if cticks is not None:
                hc.set_ticks(cticks)
            if cticklabels is not None:
                hc.ax.tick_params(labelsize=9)
                cax.set_xticklabels(cticklabels, fontsize=10)

            hc.locator = ticker.MaxNLocator(integer=True, min_n_ticks=2, nbins=5)
            hc.update_ticks()

            cax.set_ylabel(clabel, rotation=0, ha='right', va='center', fontsize=12)

    return ax
