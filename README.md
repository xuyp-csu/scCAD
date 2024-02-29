# scCAD
## Introduction
Cluster decomposition-based Anomaly Detection method (scCAD) is used to effectively identify rare cell types in scRNA-seq data.

![alt text](https://github.com/xuyp-csu/scCAD/blob/main/scCAD_overview.png)

## Getting Started
## Hardware requirements
`scCAD` package requires only a standard computer with enough RAM to support the in-memory operations.

### OS Requirements
This package is supported for *windows* and *Linux*. The package has been tested on the following systems:
+ windows: windows 10 home (1909)
+ Linux: Ubuntu 16.04.7

### Prerequisites

	Python --- 3.7.13
	h5py --- 3.7.0
	networkx --- 2.6.3
	numpy --- 1.21.5
	pandas --- 1.3.5
	python-louvain --- 0.16
	Scanpy --- 1.9.1
	scikit-learn --- 1.0.2
	scipy --- 1.7.3
	tqdm --- 4.64.0

### Installation

1. (Conda users, Conda version: 4.12.0) Create your environment and activate it:
	```
	conda create -n scCAD_env python=3.7
 	source activate scCAD_env
 	```
 	-Use `conda activate` instead of `source activate`, if required.
   
	-Building typically completes in about 30 seconds.

3. Install dependencies with pip:

	```
	pip install -r requirements.txt
	```
 	-Installation typically requires around 3 to 5 minutes, depending on network conditions.
 
## Usage
### Arguments:
	data: The expression data matrix with genes in columns and cells in rows.
	dataName: Name of scRNA-seq dataset. (default: None)
  	cellNames: Names of all cells. (default: None)
	geneNames: Names of all genes. (default: None)
	normalization: Whether the data needs to be normalized. (default: True)
	seed: Random seed. (default: 2023)
	rare_h: Rare threshold to use when doing cluster decomposition. (default: 0.01)
	merge_h: Threshold to use when doing cluster merge. (default: 0.3)
	overlap_h: Overlap threshold to identify rare sub-clusters. (default: 0.7)
	save_full: Whether the full result needs to be saved. (default: True)
	save_path:  Path to save results.
### Files:
1. 1%Jurkat.h5 -- An example scRNA-seq dataset comprising 1540 293T cells and 16 Jurkat cells.
2. scCAD.py -- The implementation of the scCAD algorithm.
3. scCAD_overview.png -- Overview of scCAD.
4. requirements.txt -- The configuration file to install the specified packages with their specified versions.
### Step-by-step description of full demo is as follows
1. Load libraries.
	```python
	import scCAD
	import numpy as np
	import pandas as pd
	import h5py
	from collections import Counter
	```
2. Load Data in current environment.
	```python
	# Data matrix should only consist of values where rows represent cells and columns represent genes.
	data_mat = h5py.File('./1%Jurkat.h5')
	data = np.array(data_mat['X']) # Cells * Genes
	labels = np.array(data_mat['Y'])
	geneNames = np.array(data_mat['gn'])
	cellNames = np.array(data_mat['cn'])
	data_mat.close()
 	data = np.vectorize(float)(data)
	labels = np.array([str(i, 'UTF-8') for i in labels])
	geneNames = np.array([str(i, 'UTF-8') for i in geneNames])
	cellNames = np.array([str(i, 'UTF-8') for i in cellNames])
	```
3. Execute scCAD on the dataset mentioned above.
	```python
 	# If gene and cell names are not provided, scCAD will generate them automatically.
	result, score, sub_clusters, degs_list = scCAD.scCAD(data=data, dataName='Jurkat', cellNames=cellNames, geneNames=geneNames, save_path='./scCAD_res/') 
 	'''
  	Returned Value :
 		result : Rare sub-clusters identified by scCAD: list.
 		score : Score of every sub_clusters: np.array[n sub-clusters].
 		sub_clusters : Assignment of sub-cluster labels for each cell: np.array[n cells].
 		degs_list : List of differentially expressed genes used for rare sub-clusters: list.
  	'''
	```
 4. View the identified results, if labels are available.
 	```python
	# If cell names are not provided, please run:
	# cellNames = [str(i) for i in range(data.shape[0])]
	for i in result:
		indices = np.where(np.isin(cellNames, i))[0]
		print(Counter(labels[indices]))
  	```

## Contact
If any questions, please do not hesitate to contact us at: 

Yunpei Xu, xu_yunpei@csu.edu.cn

Jianxin Wang, jxwang@csu.edu.cn

## How to cite?
If you use this tool, please cite the following work.

Xu Y, et al. Cluster decomposition-based Anomaly Detection for Rare Cell Identification in Single Cell Expression Data. 2023, submitted.

