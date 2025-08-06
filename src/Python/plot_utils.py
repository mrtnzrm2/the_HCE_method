import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from utils import *

def draw_dendrogram(Z : npt.NDArray, R : int, nodes, cmap_name="hls", leaf_font_size=20,  linewidth=0.5, remove_labels=False, **kwargs):
      ''' Draws a dendrogram from the linkage matrix Z.
    
      Parameters
      ----------

      Z : npt.NDArray
          Linkage matrix in Scipy format.
      R : int
          Number of clusters to visualize.
      nodes : int
          Number of nodes in the dendrogram.
      cmap_name : str, optional
          Name of the colormap to use for coloring clusters. Default is "hls".
      leaf_font_size : int, optional
          Font size for leaf labels. Default is 20.
      linewidth : float, optional
          Width of the lines in the dendrogram. Default is 0.5.
      remove_labels : bool, optional
          If True, removes labels from the leaves of the dendrogram. Default is False.
      **kwargs : dict, optional
          Additional keyword arguments to pass to the dendrogram function.
      '''

      print("Visualize Zdendrogram!!!")
      from scipy.cluster import hierarchy
      import matplotlib
      from matplotlib.colors import to_hex

      matplotlib.rcParams['lines.linewidth'] = linewidth
      fig, ax = plt.subplots(1, 1)

      if R > 1:
        partition = fast_cut_tree(Z, n_clusters=R).ravel()
        new_partition = filter_partition(partition)
        unique_clusters_id = np.sort(np.unique(new_partition))
        if -1 in unique_clusters_id:
          cm = sns.color_palette(cmap_name, len(unique_clusters_id)-1)
        else:
          cm = sns.color_palette(cmap_name, len(unique_clusters_id))
        dlf_col = "#808080"
        ##
        D_leaf_colors = {}
        for i, _ in enumerate(np.arange(nodes)):
          if new_partition[i] != -1:
            D_leaf_colors[i] = to_hex(cm[new_partition[i]])
          else: D_leaf_colors[i] = dlf_col
        ##
        link_cols = {}
        for i, i12 in enumerate(Z[:,:2].astype(int)):
          c1, c2 = (link_cols[x] if x > len(Z) else D_leaf_colors[x]
            for x in i12)
          link_cols[i+1+len(Z)] = c1 if c1 == c2 else dlf_col

        if not remove_labels:
          hierarchy.dendrogram(
            Z,
            labels=np.arange(nodes),
            color_threshold=Z[nodes - R, 2],
            link_color_func = lambda k: link_cols[k],
            leaf_rotation=90, leaf_font_size=leaf_font_size, ax=ax, **kwargs
          )
        else:
          hierarchy.dendrogram(
            Z,
            no_labels=True,
            color_threshold=Z[nodes - R, 2],
            link_color_func = lambda k: link_cols[k], leaf_rotation=90, ax=ax
          )
      else:
        deep_blue = to_hex(sns.color_palette("deep")[0])
        if not remove_labels:
          hierarchy.dendrogram(
            Z,
            labels=np.arange(nodes),
            color_threshold=np.Inf,
            link_color_func = lambda k: deep_blue,
            leaf_rotation=90, leaf_font_size=leaf_font_size, ax=ax, **kwargs
          )
        else:
          hierarchy.dendrogram(
            Z,
            no_labels=True,
            color_threshold=np.Inf,
            link_color_func = lambda k: deep_blue, leaf_rotation=90, ax=ax
          )
      plt.gca().set_ylabel("Distance")
      sns.despine()
      fig.set_figwidth(15)
      fig.set_figheight(7)
      plt.xticks(rotation=90)