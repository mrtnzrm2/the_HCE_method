# The HCE Method

Official code for the paper:  
**Hierarchical community detection via maximum entropy partitions and the renormalization group**  
by Jorge Martinez Armas

## Overview

Hierarchical community structure is a hallmark of complex systems. Many algorithms compute dendrograms to reveal this structure, but identifying meaningful levels is challenging.

We introduce the **Hierarchical Clustering Entropy (HCE)** method, combining information theory and statistical physics to optimize multilevel partitioning in dendrograms. HCE is tested on synthetic benchmarks (flat, symmetric, asymmetric hierarchies) and real-world networks, from high-school social graphs to larval zebrafish neural networks.

**Key findings:**  
- HCE reliably identifies hierarchical organization in networks with known ground truth, achieving performance comparable to state-of-the-art methods.
- It also uncovers known structural and functional organization in diverse real-world datasets.

## Repository Structure

```
.
├── src/
│   ├── Python/
│   │   └── hce_framework.py        # Main Python module for HCE
│   ├── MATLAB/
│   │   └── utils/                  # MATLAB utility functions
├── scripts/
│   ├── MATLAB/
│   │   ├── exampleHRNG.m           # HNRG (symmetric) example
│   │   └── exampleHB.m             # Hierarchical benchmark (asymmetric) example
├── notebooks/
│   └── examples.ipynb              # Python examples: Karate Club, HNRG
├── requirements.txt
└── README.md
```

## Getting Started

### Create a Python virtual environment (recommended)

It is recommended to use a virtual environment to manage dependencies.  
To create a virtual environment with Python 3.11.9, run:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
```

If you do not have Python 3.11.9 installed, download it from [python.org](https://www.python.org/downloads/release/python-3119/).

### Python

1. Clone the repository:
    ```bash
    git clone https://github.com/mrtnzrm2/the_HCE_method.git
    cd the_HCE_method
    ```
2. Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```
3. Run examples in the Jupyter notebook:
    ```bash
    jupyter notebook notebooks/examples.ipynb
    ```

### MATLAB

Requires **Statistics and Machine Learning Toolbox** (see below).

## MATLAB Prerequisites
This repository requires the **Statistics and Machine Learning Toolbox** from MATLAB.

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

## License

This project is licensed under the [BSD 3-Clause License](LICENSE).

## Contact

Contact: Jorge Martinez Armas (<jmrtnza@gmail.com>)  
For issues or support, please use the [GitHub Issues](https://github.com/mrtnzrm2/the_HCE_method.git/issues).