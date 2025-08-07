import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from utils import *

def draw_HCE(hce : dict, s=15, xshift=0, yshift=0, linewidth=1, fonttextsize=20, fontlabelsize=20, fontticklabelsize=12, color="r", f=None, xlabel=r"$K$", ylabel=r"HCE", ax=None):
  ''' Draws the HCE values for each level of the hierarchy.
  
  Parameters
  ----------
  hce : dict
      Dictionary with the HCE metric for each level of the hierarchy.
  s : int, optional
      Size of the marker for the maximum HCE value. Default is 15.
  xshift : float, optional
      Vertical shift for the text label of the maximum HCE value. Default is 0.
  yshift : float, optional
      Horizontal shift for the text label of the maximum HCE value. Default is 0.
  linewidth : float, optional
      Width of the line in the plot. Default is 1.
  color : str, optional
      Color of the marker for the maximum HCE value. Default is "r".
  f : callable, optional
      Function to apply to the x-axis values (K). If None, K is used as is.
  xlabel : str, optional
      Label for the x-axis. Default is r"$K$".
  ylabel : str, optional
      Label for the y-axis. Default is r"HCE".
  ax : matplotlib.axes.Axes, optional
      Axes object to draw the plot on. If None, uses the current axes.
      Default is None.
  '''

  K = np.array(list(hce.keys())) # K is the number of communities per level
  k_sort = np.argsort(K)
  hce_values = np.array(list(hce.values()))
  K = K[k_sort]
  hce_values = hce_values[k_sort]

  argmax_hce = np.nanargmax(hce_values)

  max_K = K[argmax_hce]
  max_hce = hce_values[argmax_hce]

  real_max_K = max_K

  if callable(f):
    K = f(K)
    f_max_K = f(max_K)
  else:
    f_max_K = max_K

  if ax is None:
    ax = plt.gca()

  ax.text(f_max_K + xshift, (max_hce + yshift), f"K = {real_max_K}", horizontalalignment='center', zorder=4, fontdict=dict(size=fonttextsize))
  ax.text(f_max_K + xshift, (max_hce + yshift), f"K = {real_max_K}", horizontalalignment='center', zorder=3, fontdict=dict(size=fonttextsize),  bbox=dict(edgecolor="k", facecolor='gray', alpha=0.1))
  ax.scatter([f_max_K], [max_hce], s=s, color=color, edgecolors='k', linewidth=linewidth, zorder=2)
  ax.plot(K, hce_values, zorder=1)

  ax.tick_params(axis='both', which='major', labelsize=fontticklabelsize)
  ax.set_ylabel(ylabel, fontdict=dict(size=fontlabelsize))
  ax.set_xlabel(xlabel, fontdict=dict(size=fontlabelsize))

  ax.minorticks_on()

  sns.despine(ax=ax)


def draw_dendrogram(Z : npt.NDArray, R : int, nodes, cmap_name="hls", leaf_font_size=20,  linewidth=0.5, remove_labels=False, ax=None, **kwargs):
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
      ax : matplotlib.axes.Axes, optional
          Axes object to draw the plot on. If None, uses the current axes.
          Default is None.
      **kwargs : dict, optional
          Additional keyword arguments to pass to the dendrogram function.
      '''

      print("Visualize Zdendrogram!!!")
      from scipy.cluster import hierarchy
      import matplotlib
      from matplotlib.colors import to_hex

      matplotlib.rcParams['lines.linewidth'] = linewidth
      if ax is None:
        ax = plt.gca()
      else:
        _, ax = plt.subplots(1, 1)

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
            color_threshold=np.inf,
            link_color_func = lambda k: deep_blue,
            leaf_rotation=90, leaf_font_size=leaf_font_size, ax=ax, **kwargs
          )
        else:
          hierarchy.dendrogram(
            Z,
            no_labels=True,
            color_threshold=np.inf,
            link_color_func = lambda k: deep_blue, leaf_rotation=90, ax=ax
          )
      plt.gca().set_ylabel("Distance")
      sns.despine()
      plt.xticks(rotation=90)

def plot_network_linkage(G : nx.Graph, H : npt.NDArray, nc : int, pos=None, node_communities=None, node_size=40, linewidth=2, ax=None):
  ''' 
  Plots a network with nodes colored by their community membership.

  Parameters
  ----------

  G : nx.Graph
      NetworkX graph object representing the network.
  H : npt.NDArray
      Hierarchical clustering linkage matrix.
  nc : int
      Number of communities to visualize.
  pos : dict, optional
      Dictionary with node positions. If None, uses Kamada-Kawai layout.
      Default is None.
  node_communities : list or np.ndarray, optional
      List or array of node community labels. If None, computes communities from H.
      Default is None.
  node_size : int, optional
      Size of the nodes in the plot. Default is 40.
  linewidth : float, optional
      Width of the edges in the plot. Default is 2.
  ax : matplotlib.axes.Axes, optional
      Axes object to draw the plot on. If None, uses the current axes.
      Default is None.
  '''

  if ax is None:
    ax = plt.gca()

  from matplotlib.colors import to_hex
  if node_communities is None:
    K = fast_cut_tree(H, n_clusters=nc)
    K = filter_partition(K)
  elif isinstance(node_communities, list) or isinstance(node_communities, np.ndarray):
    K = node_communities
  else:
    raise ValueError("node_communities is not a list.")

  unique_K = np.sort(np.unique(K))
  dft_color = to_hex((0.5, 0.5, 0.5))

  if -1 in unique_K:
    cm = list(sns.color_palette("hls", len(unique_K)-1))
    cm = [dft_color] + cm
  else:
    cm = list(sns.color_palette("hls", len(unique_K)))
    
  cm = {str(u): to_hex(c) for u, c in zip(unique_K, cm)}

  if pos is None:
    pos = nx.kamada_kawai_layout(G)
  
  nx.draw_networkx_nodes(
    G, pos=pos, node_color=[cm[str(u)] for u in K],
    node_size=node_size, ax=ax
  )

  nx.draw_networkx_edges(G, pos=pos, width=linewidth, ax=ax)
  ax.axis('off')