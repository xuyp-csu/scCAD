# scCAD
## Introduction
Cluster decomposition-based Anomaly Detection method (scCAD) is used to effectively identify rare cell types in scRNA-seq data.

![alt text](https://github.com/xuyp-csu/scCAD/blob/main/scCAD_overview.png)

## Getting Started
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

1. (Conda users) Create your environment and activate it:
	```
	conda create -n scCAD_env python=3.7
 	source activate scCAD_env
 	```  

2. Install dependencies with pip:

	```
	pip install -r requirements.txt
	```
 
## Usage
### Arguments:

data: Gene expression data matrix, gene in columns and samples in rows.

dataName: Name of scRNA-seq dataset.

cellNames: Names of all cells in the scRNA-seq dataset.

geneNames: Names of all genes in the scRNA-seq dataset.

normalization: Whether the data needs to be normalized. (default: True)

seed: Random seed.

merge_h: Threshold to use when doing cluster merge. (default: 0.3)

overlap_h: Overlap threshold to identify rare clusters. (default: 0.7)

rare_h: Rare threshold to use when doing cluster decomposition. (default: 0.01)

save_full: Whether the full result needs to be saved. (default: True)

save_path: Path to save results.

### Files:
UUOkidney.h5 -- a example scRNAseq data, UUOkidney

scCAD.py -- implementation of scCAD algorithm

scCAD_overview.png -- CellBRF workflow

### Step-by-step description of full demo is as follows
1. <h4>Load libraries </h4>

	```python
	import scCAD
	import numpy as np
	import pandas as pd
	import h5py
	from collections import Counter
	```
2. <h4>Load Data in current environment.</h4>
	```python
	# Data matrix should only consist of values where rows represent cells and columns represent genes.
	data_mat = h5py.File('./1%Jurkat.h5')
	data = np.array(data_mat['X']) # Cells * Genes
	labels = np.array(data_mat['Y'])
	geneName = np.array(data_mat['gn'])
	cellName = np.array(data_mat['cn'])
	data_mat.close()
	labels = np.array([str(i, 'UTF-8') for i in y])
	geneName = np.array([str(i, 'UTF-8') for i in geneName])
	cellName = np.array([str(i, 'UTF-8') for i in cellName])
	```
3. <h4>Execute scCAD on the dataset mentioned above.</h4>
	```python
 	# If gene and cell names are not provided, scCAD will generate them automatically.
	result, score, sub_clusters, degs_list = scCAD.scCAD(data=data, dataName='Jurkat', cellNames=cellName, geneNames=geneName, save_path='./scCAD_res/') 
 	'''
  	Returned Value :
		result : Rare subclusters identified by scCAD: list.
     	    	score : Score of every sub_clusters: np.array[n sub_clusters].
	    	sub_clusters : Sub-cluster label assignment for each cell: np.array[n cells].
      		degs_list : List of differentially expressed genes used for rare sub-clusters: list.
  	'''
	```
 4. <h4>View the identified results, if labels are available.</h4>
 	```python
  	
  	```

## Contact
If any questions, please do not hesitate to contact us at: 

Yunpei Xu, xu_yunpei@csu.edu.cn

Jianxin Wang, jxwang@csu.edu.cn

## How to cite?
If you use this tool, please cite the following work.

Xu Y, et al. Cluster decomposition-based Anomaly Detection for Rare Cell Identification in Single Cell Expression Data. 2023, submitted.

