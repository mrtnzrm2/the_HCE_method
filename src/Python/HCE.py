import numpy as np
import numpy.typing as npt
from utils import *

def HCE(H : npt.NDArray, N : int, use_tqdm=False):
    '''
    Computes the hierarchical clustering of a dendrogram represented by a linkage matrix.
    Parameters
    ----------
    H : npt.NDArray
        Hierarchical clustering linkage matrix.
    N : int
        Number of nodes in the hierarchy.
    use_tqdm : bool, optional
        If True, uses tqdm for progress bar.
    Returns
    -------
    hce : dict
        Dictionary with the HCE metric for each level of the hierarchy.
    '''

    if H.shape[0] != N - 1:
        raise ValueError(f"Expected {N-1} rows in H, got {H.shape[0]}.")
    if H.shape[1] != 3:
        raise ValueError(f"Expected 3 columns in H, got {H.shape[1]}.")
    
    from tqdm import tqdm
    
    # Initialize the tree structure
    # Each node will store its size
    # This is a dictionary where keys are node indices and values are dictionaries with size
    # This allows for dynamic size updates as we merge nodes
    # T[(i)] = {"size": size_of_node_i}
    # Size is initialized to 0 for all nodes for node removal of the HCE process

    T = {(i) : {"size" : np.array(0)} for i in np.arange(N)}

    hce = {}

    if use_tqdm:
       iterator = tqdm(np.arange(N-1))
    else:
       iterator = np.arange(N-1)

    for i in iterator:
        nx, ny, h = int(H[i, 0]), int(H[i, 1]), H[i, 2]
        
        T[(N + i)] = {"size" : T[(nx)]["size"] + T[(ny)]["size"] + 1}   # Update size of the new node

        del T[(nx)]
        del T[(ny)]

        s = 0 
        for key in T.keys():
            if T[key]["size"] == 0: continue
            s += -(T[key]["size"] / (i+1)) * np.log(T[key]["size"]/ (i+1))  # Compute entropy
        
        if i < N - 2:
          if h != H[i+1, 2]:
           hce[(N-(i+1))] = s * (i+1) / (N-1) # apply community size normalization
        else:
          hce[(N-(i+1))] = s * (i+1) / (N-1) # apply community size normalization

    return hce


def rHCE(H : npt.NDArray, N : int, rN : int, use_tqdm=False):
    '''
    Computes the renormalized hierarchical clustering of a dendrogram represented by a linkage matrix.
    Parameters
    ----------
    H : npt.NDArray
        Hierarchical clustering linkage matrix.
    N : int
        Number of nodes in the hierarchy.
    rN : int
        Number of communities at the renormalization level..
    use_tqdm : bool, optional
        If True, uses tqdm for progress bar.
    Returns
    -------
    rhce : dict
        Dictionary with the renormalized HCE metric for each level of the hierarchy.
    '''

    if H.shape[0] != N - 1:
        raise ValueError(f"Expected {N-1} rows in H, got {H.shape[0]}.")
    if H.shape[1] != 3:
        raise ValueError(f"Expected 3 columns in H, got {H.shape[1]}.")
    if rN < 1 or rN >= N:
        raise ValueError(f"rN must be in the range [1, {N-1}]. Got {rN}.")
    
    from tqdm import tqdm

    # Initialize the tree structure
    # Each node will store its size
    # This is a dictionary where keys are node indices and values are dictionaries with size
    # This allows for dynamic size updates as we merge nodes
    # T[(i)] = {"size": size_of_node_i}
    # Size is initialized to 0 for all nodes for node removal of the HCE process

    T = {(i) : {"size" : np.array(0)} for i in np.arange(N)}

    rhce = {}

    # Update the tree structure for the first N - rN nodes
    # These nodes will be removed in the renormalization process
    # We initialize their size to 0, as they will not contribute to the renormalization

    for i in np.arange(N - rN):
        nx, ny = int(H[i, 0]), int(H[i, 1])
        
        T[(N + i)] = {"size" : np.array(0)}

        del T[(nx)]
        del T[(ny)]

        rhce[(N-(i+1))] = 0

    if use_tqdm:
       iterator = tqdm(enumerate(np.arange(N - rN, N - 1)))
    else:
       iterator = enumerate(np.arange(N - rN, N - 1))

    for j, i in iterator:
        nx, ny, h = int(H[i, 0]), int(H[i, 1]), H[i, 2]
        
        T[(N + i)] = {"size" : T[(nx)]["size"] + T[(ny)]["size"] + 1}

        del T[(nx)]
        del T[(ny)]

        s = 0 
        for key in T.keys():
            if T[key]["size"] == 0: continue
            s += -(T[key]["size"] / (j+1)) * np.log(T[key]["size"]/ (j+1))  # Compute entropy
        
        if i < N - 2:
          if h != H[i+1, 2]:
            rhce[(N-(i+1))] = s * (j+1) / (rN-1)    # apply community size normalization
        else:
            rhce[(N-(i+1))] = s * (j+1) / (rN-1)    # apply community size normalization

    return rhce