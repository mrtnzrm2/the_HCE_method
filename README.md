# The HCE method

Official code for the paper:

**Hierarchical community detection via maximum entropy partitions and the renormalization group**

by Jorge Martinez Armas

## Overview
- Complex systems often feature a community and hierarchical organization. Several algorithms exist to compute dendrograms encoding the hierarchical community structure of a system. However, it is challenging to identify suitable dendrogram levels with meaningful clustering information.
- Our contribution is to develop the *Hierarchical Clustering Entropy* (HCE) framework that combines elements from information theory :floppy_disk: and physics :rocket: to optimize the discovery of multilevel partitions of a dendrogram to get a more comprenhensive of the interation of elemets of a system at multiple resolutions.
- The framework is tested in several benchmarks exhibiting flat, symmetric, and asymmetric hierarchical structure, as well as, real-world networks with hundreds (high-school :school:) to thoursands (larval zebrafish :fish:) nodes.
- Our key findings are that HCE (1) identifies with state-of-the-art accuracy the hierarhical organization of networks with ground truth, and (2) known structural and functional organizations in reald-world networks.

## Repository structure

<pre lang="markdown"><code>
.
├── src/
│   ├── Python/
│       ├── hce_framework.py    # Main
│   ├── MATLAB/
│       ├── utils/
├── scripts/
│   ├── MATLAB/
│       ├── exampleHRNG.m       # hierarchical nested random graph (HNRG; symmetric) example in MATLAB
│       ├── exampleHB.m         # hierarhical benchmark (HB; asymmetric) example in MATLAB
├── notebooks/
│   ├── examples.ipynb          # Karate Club and HNRG model examples in Python
├── requirements.txt
└── README.md
</code></pre>


## MATLAB Prerequisites
This repository requires the **Statistics and Machine Learning Toolbox** from MATLAB.

### Why do I need this?
It is needed to run the HNRG model.

### How to install
You can install the toolbox using one of the following methods:

- **MATLAB Add-On Explorer:**  
  1. Open MATLAB.
  2. Go to the "Home" tab.
  3. Click "Add-Ons" > "Get Add-Ons".
  4. Search for **Statistics and Machine Learning Toolbox** and install.

- **MATLAB Command Line:**  
  Run the following command in MATLAB:
  ```matlab
  matlab.addons.install('Statistics and Machine Learning Toolbox')
  ```

> **Note:** Ensure you are using MATLAB R2021a or newer for full compatibility.