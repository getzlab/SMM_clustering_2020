import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from typing import Union
import pandas as pd
import signatureanalyzer as sa

# Cluster Color map
COLORMAP={
    0:'#0073C2FF',
    1:'#EFC000FF',
    2:'#868686FF',
    3:'#CD534CFF',
    4:'#7AA6DCFF',
    5:'#003C67FF',
}

COLORMAP2={
    "C1":'#0073C2FF',
    "C2":'#EFC000FF',
    "C3":'#868686FF',
    "C4":'#CD534CFF',
    "C5":'#7AA6DCFF',
    "C6":'#003C67FF'
}

# Plotting funcs
def plot_scatter(df: pd.DataFrame, group: str, ax: plt.Axes, pal=None):
    """
    Generic scatterplot.
    ------------------------------
    Inputs:
        * df: dataframe wtih columns as axes
        * group: column in dataframe that defines group colors
        * ax: matplotlib axes
        * pal: dictionary mapping groups to colors

    Output:
        * none
    """
    groups = list(set(df[group]))

    if pal is None:
        groups_cdict = {groups[i]:x for i,x, in enumerate(sns.color_palette("hls", len(set(df[group]))))}
    else:
        groups_cdict = pal

    _ = ax.scatter(
        df['PC1'],
        df['PC2'],
        alpha=0.8,
        c=df[group].apply(lambda x: groups_cdict[x]),
        label=None,
        s=70,
        edgecolor='black',
        linewidth=0.5
    )

    ax.set_xlabel("PC1", fontsize=14)
    ax.set_ylabel("PC2", fontsize=14)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    for k,v in groups_cdict.items():
        try:
            m = ax.scatter(-1e4, -1e4, alpha=.8, c=np.array(v)[np.newaxis,:], label=k)
        except:
            m = ax.scatter(-1e4, -1e4, alpha=.8, c=v, label=k)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.legend(loc='lower left')
    ax.set_title(group, fontsize=14, fontweight='bold', x=0.05, y=0.025, ha='left')

def plot_metric_hm(X, sample_n, k_n, title: str = None, figsize: tuple = (8,6), ax: plt.Axes = None):
    """
    Metrics heatmap for downsampling analysis.
    ------------------------------
    Inputs:
        * X
        * sample_n
        * k_n
        * title:
        * figsize: tuple of
        * ax: plt.Axes

    Outputs:
        * none
    """
    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        pd.DataFrame(X.mean(2), index=["n={}".format(str(x)) for x in sample_n], columns=[str(x) for x in k_n]),
        annot=True,
        cmap='coolwarm',
        ax=ax,
        linewidth=0.1,
        linecolor='black',
        fmt='.2f'
    )

    [ax.spines[sp].set_visible(True) for sp in ['top','right','left','bottom']]

    ax.set_yticklabels(ax.get_yticklabels(), rotation = 0, fontsize=12, fontstyle='italic')
    ax.set_xlabel("K Factors", fontsize=16)

    if title is not None:
        ax.set_title(title, fontsize=18)

def plot_dist_per_metric(X, k, sample_n, k_n, figsize=(8,6), ax=None, title=None, s=10):
    """
    Plot distribution per for a given K.
    """
    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)

    X = pd.DataFrame(X[:,k-min(k_n),:], index=["n={}".format(str(x)) for x in sample_n]).T

    sns.stripplot(
        data=X,
        s=s,
        ax=ax,
        linewidth=.25,
        alpha=0.25,
    )

    sns.violinplot(
        data=X,
        ax=ax,
        linewidth=1,
        alpha=0.6,
        color='White'
    )

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=14)
    ax.set_title("K = {}".format(k), fontsize=14)

def plot_dist_per_metric_grid(X, sample_n, k_n, y_label=None):
    """
    Plot distribution grid for a given metric.
    """

    fig,axes = plt.subplots(3,3,figsize=(12,12), sharex=True)

    c=0
    for i in range(3):
        for j in range(3):
            plot_dist_per_metric(X, k_n[c], sample_n, k_n, ax=axes[i,j], s=8)
            c+=1
            if j==0: axes[i,j].set_ylabel(y_label, fontsize=14)

    plt.tight_layout()

def plot_marker_heatmap(
    X: pd.DataFrame,
    signatures: pd.DataFrame,
    order_series: pd.Series,
    signatures_idx: str = 'max_id',
    subset_genes: Union[pd.Series,None] = None,
    diff: float = 0.5,
    max_norm: float = 0.5,
    figsize: tuple = (16,12),
    cmap: str ="YlGnBu",
    vmax: float = None,
    vmin: float = None,
    cohort_s: Union[pd.Series,None] = None
    ):
    """
    Plot marker map.
    -----------------------------
    Args:
        * X: pd.DataFrame of input sample x feature matrix
        * signatures: pd.DataFrame signatures output;
            this bundles information about the weightings of each feature (ex. gene) and
            what signature they map to
        * order_series: series of samples mapping to subgroups
            index: X.index
            values: subgrouping
        * signatures_idx: string for signature grouping
        * subset_series: a pd.Series with the index as the gene name or ID that
            matches the marker matrix & has a "Subgroup" column for labeling
        * diff: difference of loading for called signature vs. rest
        * max_norm: strength of loading for called signature
        * figsize: size of figure
        * cmap: colormap for plot
        * display_y: whether or not to display feature names
        * vmax: colorbar max
        * vmin: colorbar min
        * cohort_s: cohort_series dataframe (added on top of plot)
        * y_hm_label: label of y-axis on heatmap (ex. Genes, Protein LFC, etc.)
        * cbar_hm_label: label of heatmap colorbar
    Returns:
        * plt.Figure
    """
    from scipy.cluster import hierarchy
    import scipy.cluster.hierarchy as shc
    from sklearn.cluster import AgglomerativeClustering

    import signatureanalyzer as sa

    # Remove signatures with no marker genes associated
    order_series = order_series[order_series.isin(set(signatures[signatures_idx].astype(int)))]

    # Filter X matrix
    sample_markers = X.loc[signatures.index, order_series.sort_values().index]

    # Set horizontal lines
    hz_lines = np.unique(sample_markers.join(signatures).loc[:,signatures_idx].values, return_index=True)[1]

    fig, ax = plt.subplots(figsize=figsize)

    x0 = ax.get_position().x0
    x1 = ax.get_position().x1
    y0 = ax.get_position().y0
    y1 = ax.get_position().y1
    buf = y1*0.01

    sns.heatmap(sample_markers, ax=ax, cmap=cmap, rasterized=True, vmax=vmax, vmin=vmin, cbar=False)
    v,c = np.unique(order_series, return_counts=True)

    # plot horizontal lines
    _c = np.cumsum(c)
    _ci = np.roll(_c,2)
    _ci[0] = 0
    _ci[1] = 0
    ax.hlines(hz_lines, _ci, _c, rasterized=True)

    # plot vertical lines
    _h = list(hz_lines)
    _h.append(ax.get_ylim()[0])
    ax.vlines(np.cumsum(c)[:-1], _h[:-2], _h[2:], rasterized=True)
    ax.vlines(np.cumsum(c)[:-1], *ax.get_ylim(), alpha=0.4, rasterized=True, linewidth=1)

    # Set yticks
    ax.yaxis.tick_right()
    ax.set_yticks(np.arange(sample_markers.index.values.shape[0])+0.5)
    ax.set_yticklabels(sample_markers.index.values, fontsize=7.5, rasterized=True, rotation=0, va="center")

    # --------------cluster annot-------------------
    clust_ax = fig.add_axes([x0, y1+buf, x1*.861, 2*buf])

    clust_ax.set_xticks([])
    clust_ax.set_yticks([])

    colors_conversion, meta_colormap = sa.pl.series_to_colors(order_series.loc[sample_markers.columns])
    meta_colormap_inv = dict([[v,k] for k,v in meta_colormap.items()])
    meta_colormap_inv = {(k[0],k[1],k[2]):v for k,v in meta_colormap_inv.items()}

    mat,cmap = sa.pl.color_list_to_matrix_and_cmap(colors_conversion)

    sns.heatmap(
        mat,
        cmap=cmap,
        ax=clust_ax,
        yticklabels=False,
        xticklabels=False,
        cbar=False
    )

    [spine.set_visible(True) for _, spine in clust_ax.spines.items()]

    clust_ax.yaxis.set_label_position("right")
    clust_ax.set_ylabel("Consensus NMF", rotation=0, va='center', ha='left')
    # --------------cluster annot-------------------


    # --------------sample annot-------------------
    if cohort_s is not None:
        cdict = {'Low': 'Green', 'Intermediate':'Yellow', 'High': 'Red'}
        order_dict = {'Green': 0, 'Yellow': 1, 'Red': 2}

        # Get ordering and samples
        cohort_s = cohort_s.loc[sample_markers.columns]

        # Create axis
        cs_ax = fig.add_axes([x0, y1+4*buf, x1*.861, 2*buf])
        cs_ax.set_xticks([])
        cs_ax.set_yticks([])

        cbar_cs_ax = fig.add_axes([x0, y1+7*buf, x1*.25, 2*buf])

        colors_conversion, meta_colormap = sa.pl.series_to_colors(cohort_s, cdict=cdict)
        meta_colormap_inv = dict([[v,k] for k,v in meta_colormap.items()])

        mat,cmap = sa.pl.color_list_to_matrix_and_cmap(colors_conversion, order_dict=order_dict)

        sns.heatmap(
            mat,
            cmap=cmap,
            ax=cs_ax,
            yticklabels=False,
            xticklabels=False,
            cbar=True,
            cbar_ax=cbar_cs_ax,
            cbar_kws={"orientation": "horizontal"}
        )

        cb_ticks = [float(t.get_text().replace('−','-')) for t in cbar_cs_ax.get_yticklabels()]

        color_value_mapping = dict([[v,k] for k,v in order_dict.items()])

        cbar_cs_ax.get_xaxis().set_ticks([])

        n_labels = len(list(color_value_mapping.keys()))

        # FIX THIS
        vals = [x * ((n_labels)/(n_labels+1)) + 0.5 * ((n_labels)/(n_labels+1)) for x in list(color_value_mapping.keys())]
        #cbar_cs_ax.get_xaxis().set_ticks(vals)

        cbar_cs_ax.get_xaxis().set_ticks([0.375, 1, 1.675])
        cbar_cs_ax.get_xaxis().set_ticklabels(list(cdict.keys()))
        cbar_cs_ax.xaxis.set_ticks_position('top')

        cbar_cs_ax.set_frame_on(True)
        [spine.set_visible(True) for _, spine in cs_ax.spines.items()]

        cs_ax.yaxis.set_label_position("right")
        cs_ax.set_ylabel("Risk", rotation=0, va='center', ha='left')

    # --------------sample annot-------------------

    # --------------pval barplot-------------------
    p_ax = fig.add_axes([x1+12*buf, y0, 10*buf, y1-y0])
    p_ax.set_yticks([])

    log10_pval_adj = -np.log10(signatures.loc[sample_markers.index]['pval_adj'])

    p_ax.barh(np.arange(signatures.shape[0]), log10_pval_adj[::-1], edgecolor='black', linewidth=1, color='darkblue')
    plt.margins(y=0)
    p_ax.axvline(1, linewidth=1, color='red')

    p_ax.spines['top'].set_visible(False)
    p_ax.spines['right'].set_visible(False)
    p_ax.set_xticks([0,5,10,15,20])
    p_ax.set_xlabel("$-log_{10}$ (adj. p-val)")
    # --------------pval barplot-------------------

    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('')

    [spine.set_visible(True) for _, spine in ax.spines.items()]

    # Set xticks
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_xticks(np.cumsum(c)-c/2)
    ax.set_xticklabels(v, rotation=360, fontsize=14)
    ax.tick_params(axis='x', which=u'both',length=0)

    return fig

def plot_consensus_matrix(
    cmatrix: pd.DataFrame,
    metric: str = 'euclidean',
    method: str = 'ward',
    n_clusters: int = 10,
    color_thresh_scale: float = 0.3,
    figsize: tuple = (8,8),
    p: int = 30,
    metas: Union[list, None] = None,
    vmax: Union[float, None] = None,
    vmin: Union[float, None] = None,
    cbar_label: str = 'ARD-NMF \nMembership',
    cmap: Union[str, None] = None,
    plot_cluster_lines: bool = False
    ):
    """
    Plot consensus matrix.
    -----------------------
    Args:
        * cmatrix: consensus matrix. This may be generated by calling:
            df, assign_p = consensus_cluster_ardnmf(filepath)
        * metric: distance metric
        * method: method of clustering
        * n_clusters: number of clusters for agglomerative clustering
        * color_thresh_scale: asthetic scale for coloring of dendrogram
        * figsize: figsize
        * p: parameter for dendrogram
        * meta: list of pd.Series that includes a variable of interest to plot
            to left of plot; must be categorical in nature
    Returns:
        * fig
    """
    from matplotlib.pyplot import cm
    import matplotlib as mpl
    from scipy.cluster import hierarchy
    import scipy.cluster.hierarchy as shc
    from sklearn.cluster import AgglomerativeClustering
    # -------------
    # Heatmap
    # -------------
    fig,ax = plt.subplots(figsize=figsize)
    cbar_ax = fig.add_axes([ax.get_position().x1 + ax.get_position().x1*0.1, ax.get_position().y0, .025, .1])

    # Compute initial linkage to grab ordering
    d_linkage = shc.linkage(cmatrix, metric=metric, method=method)
    dres = shc.dendrogram(d_linkage, p=p, no_plot=True)
    dgram_idx = list(map(int, dres['ivl']))

    # Create heatmap
    if vmax is None:
        cbar_top_lim = np.max(cmatrix.values)
    else:
        cbar_top_lim = vmax

    if vmin is None:
        cbar_bottom_lim = 0
    else:
        cbar_bottom_lim = vmin

    # Create heatmap
    sns.heatmap(
        cmatrix.iloc[dgram_idx,dgram_idx].values,
        ax=ax,
        square=True,
        cbar_ax=cbar_ax,
        cbar_kws = {'ticks':[cbar_bottom_lim, cbar_top_lim]},
        rasterized=True,
        vmax=vmax,
        vmin=vmin,
        cmap=cmap
    )

    cbar_ax.set_ylabel(cbar_label, fontsize=10,rotation=90)
    ax.set_xticks([])
    ax.set_yticks([])

    x0 = ax.get_position().x0
    x1 = ax.get_position().x1
    y0 = ax.get_position().y0
    y1 = ax.get_position().y1

    buf = y1*0.015

    # -------------
    # Clustering
    # -------------
    cluster = AgglomerativeClustering(
        n_clusters=n_clusters,
        affinity=metric,
        linkage=method
    )

    clusters = cluster.fit_predict(cmatrix.iloc[dgram_idx,dgram_idx])
    cluster_color_list, _ = sa.pl.series_to_colors(pd.Series(clusters), cdict=COLORMAP)

    # -------------
    # Dendrogram
    # -------------
    cmap = cm.rainbow(np.linspace(0, 1, 10))
    hierarchy.set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in cmap])

    dax = fig.add_axes([x0, y1+buf, x1-x0, 0.15])

    dres = shc.dendrogram(
        d_linkage,
        p=p,
        ax=dax,
        above_threshold_color="grey",
        color_threshold=color_thresh_scale*max(d_linkage[:,2])
    )

    dax.set_xticks([])
    dax.set_yticks([])
    [dax.spines[x].set_visible(False) for x in ['top','right','bottom','left']]

    # -------------
    # Metadata Axes
    # -------------
    if plot_cluster_lines:
        hz_lines = np.sort(np.unique(pd.Series(clusters), return_index=True)[1])
        v,c = np.unique(clusters, return_counts=True)

        _c = hz_lines
        _c = np.roll(hz_lines, 1)
        _c[0] = 0
        _c[1] = 0

        _ci = hz_lines[1:]
        _ci = np.append(_ci, clusters.shape[0])

        for idx, hz in enumerate(hz_lines):
            ax.hlines(hz, _c[idx], _ci[idx], rasterized=True)
            ax.vlines(hz, _c[idx], _ci[idx], rasterized=True)

    # Add axes
    # Plots agglomerative clustering results
    if metas is None:
        lax = fig.add_axes([x0-3*buf, y0, 2*buf, y1-y0])
        mat, cmap = sa.pl.color_list_to_matrix_and_cmap(cluster_color_list)
        sns.heatmap(mat.T, cmap=cmap, ax=lax, xticklabels=False, yticklabels=False, cbar=False, rasterized=True)

        uniq, idx, num_vals = np.unique(clusters.T, return_index=True, return_counts=True)
        y_locs = idx + num_vals / 2

        for idx,u in enumerate(uniq):
            lax.text(x0-50*buf, y_locs[idx], u, ha='center')

        for idx,u in enumerate(uniq):
            ax.text(
                mat.shape[1]+0.01*mat.shape[1],
                y_locs[idx],
                "n={}".format(num_vals[idx]),
                ha='left',
                fontsize=14
            )

        for _, spine in lax.spines.items():
            spine.set_visible(True)

        lax.set_xlabel("Consensus", rotation=90)

    else:
        for idx,meta in enumerate(metas):
            new_ax = [x0-(idx+3)*buf-(idx*2)*buf, y0, 2*buf, y1-y0]
            lax = fig.add_axes(new_ax)

            if isinstance(meta, str) and meta=='aggr':
                mat, cmap = sa.pl.color_list_to_matrix_and_cmap(cluster_color_list)
                sns.heatmap(mat.T, cmap=cmap, ax=lax, xticklabels=False, yticklabels=False, cbar=False, rasterized=True)

                uniq, idx, num_vals = np.unique(clusters.T, return_index=True, return_counts=True)
                y_locs = idx + num_vals / 2

                for idx,u in enumerate(uniq):
                    lax.text(0.5, y_locs[idx], "C{}".format(u+1), ha='center', color='white')

                for idx,u in enumerate(uniq):
                    ax.text(
                        mat.shape[1]+0.01*mat.shape[1],
                        y_locs[idx],
                        "n={}".format(num_vals[idx]),
                        ha='left',
                        fontsize=14
                    )

                #lax.set_xlabel("Consensus", rotation=90)

            else:
                meta = meta.loc[cmatrix.index[dgram_idx]].fillna(0).astype(int)
                cdict={1:'purple',0:'white'}

                cluster_color_list, _ = sa.pl.series_to_colors(meta, cdict=cdict)
                mat,cmap = sa.pl.color_list_to_matrix_and_cmap(cluster_color_list)
                sns.heatmap(mat.T, cmap=cmap, ax=lax, yticklabels=False, xticklabels=False, cbar=False)
                lax.set_xlabel(meta.name, rotation=90)

            for _, spine in lax.spines.items():
                spine.set_visible(True)

    rs = pd.DataFrame(clusters, index=cmatrix.index[dgram_idx]).rename(columns={0:'clusters'})

    for _, spine in ax.spines.items():
        spine.set_visible(True)

    ax.set_xlabel("Samples", fontsize=14)

    return fig, rs

def plot_marker_heatmap_fig1(
    X: pd.DataFrame,
    signatures: pd.DataFrame,
    order_series: pd.Series,
    signatures_idx: str = 'max_id',
    figsize: tuple = (16,13),
    vmax: float = None,
    vmin: float = None,
    metas: Union[pd.Series,None] = None,
    order_x: Union[pd.Series,None] = None,
    ):
    """
    Plot marker map.
    -----------------------------
    Args:
        * X: pd.DataFrame of input sample x feature matrix
        * signatures: pd.DataFrame signatures output;
            this bundles information about the weightings of each feature (ex. gene) and
            what signature they map to
        * order_series: series of samples mapping to subgroups
            index: X.index
            values: subgrouping
        * signatures_idx: string for signature grouping
        * subset_series: a pd.Series with the index as the gene name or ID that
            matches the marker matrix & has a "Subgroup" column for labeling
        * diff: difference of loading for called signature vs. rest
        * max_norm: strength of loading for called signature
        * figsize: size of figure
        * cmap: colormap for plot
        * display_y: whether or not to display feature names
        * vmax: colorbar max
        * vmin: colorbar min
        * cohort_s: cohort_series dataframe (added on top of plot)
        * y_hm_label: label of y-axis on heatmap (ex. Genes, Protein LFC, etc.)
        * cbar_hm_label: label of heatmap colorbar
    Returns:
        * plt.Figure
    """
    from matplotlib.colors import ListedColormap
    from scipy.cluster import hierarchy
    import scipy.cluster.hierarchy as shc
    from sklearn.cluster import AgglomerativeClustering

    cmap = ListedColormap(['white','darkblue'])

    # Remove signatures with no marker genes associated
    order_series = order_series[order_series.isin(set(signatures[signatures_idx]))]

    # Filter X matrix
    if order_x is None:
        order_series = order_series.sort_values()
    else:
        order_series = order_series.loc[order_x]

    sample_markers = X.loc[signatures.index, order_series.sort_values().index]

    # Set horizontal lines
    hz_lines = np.unique(sample_markers.join(signatures).loc[:,signatures_idx].values, return_index=True)[1]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    x0 = ax.get_position().x0
    x1 = ax.get_position().x1
    y0 = ax.get_position().y0
    y1 = ax.get_position().y1
    buf = y1*0.01

    sns.heatmap(sample_markers, ax=ax, cmap=cmap, rasterized=True, vmax=vmax, vmin=vmin, cbar=False)
    v,c = np.unique(order_series, return_counts=True)

    # plot horizontal lines
    _c = np.cumsum(c)
    _ci = np.roll(_c,2)
    _ci[0] = 0
    _ci[1] = 0
    ax.hlines(hz_lines, _ci, _c, rasterized=True, color='k')

    # plot vertical lines
    _h = list(hz_lines)
    _h.append(ax.get_ylim()[0])
    ax.vlines(np.cumsum(c)[:-1], _h[:-2], _h[2:], rasterized=True , color='k')
    ax.vlines(np.cumsum(c)[:-1], *ax.get_ylim(), alpha=0.8, rasterized=True, linewidth=1 , color='lightgrey')

    # Set yticks
    ax.yaxis.tick_right()
    ax.set_yticks(np.arange(sample_markers.index.values.shape[0])+0.5)
    ax.set_yticklabels(sample_markers.index.values, fontsize=7.5, rasterized=True, rotation=0, va="center")

    # --------------cluster annot-------------------
    for idx,meta in enumerate(metas):
        if meta.unique().shape[0]==2 and meta.dtype=='int64':
            meta = meta.astype(bool)
            cdict={True:'purple',False:'white'}
        elif meta.name =='Consensus':
            cdict=COLORMAP2
        else:
            cdict=None

        new_ax = [x0, y1+buf*(idx*3+1), x1*.861, 2*buf]
        clust_ax = fig.add_axes(new_ax)

        clust_ax.set_xticks([])
        clust_ax.set_yticks([])

        colors_conversion, meta_colormap = sa.pl.series_to_colors(meta.loc[sample_markers.columns],cdict=cdict)
        meta_colormap_inv = dict([[v,k] for k,v in meta_colormap.items()])
        meta_colormap_inv = {(k[0],k[1],k[2]):v for k,v in meta_colormap_inv.items()}

        mat,cmap = sa.pl.color_list_to_matrix_and_cmap(colors_conversion)

        sns.heatmap(
            mat,
            cmap=cmap,
            ax=clust_ax,
            yticklabels=False,
            xticklabels=False,
            cbar=False
        )

        [spine.set_visible(True) for _, spine in clust_ax.spines.items()]

        clust_ax.yaxis.set_label_position("right")
        clust_ax.set_ylabel(meta.name, rotation=0, va='center', ha='left')
    # --------------cluster annot-------------------

    # --------------pval barplot-------------------
    p_ax = fig.add_axes([x1+12*buf, y0, 10*buf, y1-y0])
    p_ax.set_yticks([])

    log10_pval_adj = -np.log10(signatures.loc[sample_markers.index]['pval_adj'])

    p_ax.barh(np.arange(signatures.shape[0]), log10_pval_adj[::-1], edgecolor='black', linewidth=1, color='purple')
    plt.margins(y=0)
    p_ax.axvline(1, linewidth=1, color='red')

    p_ax.spines['top'].set_visible(False)
    p_ax.spines['right'].set_visible(False)
    p_ax.set_xticks([0,5,10,15,20])
    p_ax.set_xlabel("$-log_{10}$ (adj. p-val)")
    # --------------pval barplot-------------------

    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('')

    [spine.set_visible(True) for _, spine in ax.spines.items()]

    # Set xticks
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_xticks(np.cumsum(c)-c/2)
    v = ["{}\n(n={})".format(x,y) for x,y in zip(*np.unique(order_series, return_counts=True))]
    ax.set_xticklabels(v, rotation=360, fontsize=12)
    ax.tick_params(axis='x', which=u'both',length=0)

    return fig
