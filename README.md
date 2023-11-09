# scCAD
## Introduction
Cluster decomposition-based Anomaly Detection method (scCAD) is used to effectively identify rare cell types in scRNA-seq data.

![alt text](https://github.com/xuyp-csu/scCAD/blob/main/scCAD_overview.png)

## Installation
### Required python modules

Python --- 3.7.13

Scanpy --- 1.9.1

scipy --- 1.7.3

scikit-learn --- 1.0.2

python-louvain --- 0.16

networkx --- 2.6.3

tqdm --- 4.64.0

### Conda users 

1. Create your environment and activate it:
	```
	conda create -n scCAD_env python=3.7
 	conda activate scCAD_env
 	```  

2. Install dependencies with pip:

	```
	pip install -r requirements.txt
	```
### 


 
## Usage


## Arguments:

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

## Files:
UUOkidney.h5 -- a example scRNAseq data, UUOkidney

scCAD.py -- implementation of scCAD algorithm

scCAD_overview.png -- CellBRF workflow

## Demo:
```
import scDOD
import os
import numpy as np
import pandas as pd
import h5py
from collections import Counter

dir = './h5data'
dataName = 'UUOkidney'
data_mat = h5py.File(os.path.join(dir, dataName + '.h5'))
X = np.array(data_mat['X'])
y = np.array(data_mat['Y'])
geneName = np.array(data_mat['gn'])
cellName = np.array(data_mat['cn'])
data_mat.close()
_, score, sub_clusters, _ = scDOD.scDOD(data=X, dataName=dataName, cellNames=cellName, geneNames=geneName, save_path='./scDOD_res/')
for i in range(len(np.unique(sub_clusters))):
    if score[i] >= 0.7:
        print(Counter(y[np.where(sub_clusters==i)[0]]))    
print(pd.DataFrame(y).value_counts())
```

# Contact
If any questions, please do not hesitate to contact us at: 

Yunpei Xu, xu_yunpei@csu.edu.cn

Jianxin Wang, jxwang@csu.edu.cn

# How to cite?
If you use this tool, please cite the following work.

Xu Y, et al. Cluster decomposition-based Anomaly Detection for Rare Cell Identification in Single Cell Expression Data. 2023, submitted.

