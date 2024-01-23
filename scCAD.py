"""
MIT License

Copyright (c) [year] [fullname]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import scanpy as sc
import numpy as np
import networkx as nx
from community import community_louvain
from sklearn.ensemble import IsolationForest, RandomForestClassifier
from sklearn.decomposition import PCA
from collections import Counter
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm
import warnings
import time
# from utils_functions import *

def normalize(adata, filter_min_counts=True, size_factors=True, normalize_input=False, logtrans_input=True):
    if filter_min_counts:
        sc.pp.filter_genes(adata, min_cells=3)
#         sc.pp.filter_cells(adata, min_genes=200)
    if size_factors or normalize_input or logtrans_input:
        adata.raw = adata.copy()
    if size_factors:
        sc.pp.normalize_per_cell(adata)
        adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
    else:
        adata.obs['size_factors'] = 1.0
    if logtrans_input:
        sc.pp.log1p(adata)
    if normalize_input:
        sc.pp.scale(adata)
    return adata

def fast_clustering(data, k=15, seed=2023):
    # scaler = StandardScaler()
    # data = scaler.fit_transform(data)
    nps = min(40, data.shape[0])
    if data.shape[0] >= nps:
        pca = PCA(n_components=nps, random_state=seed)
        X_pca = pca.fit_transform(data)
    else:
        X_pca = data

    from sklearn.neighbors import NearestNeighbors
    nn_res = NearestNeighbors(n_neighbors=(k+1), algorithm='auto', n_jobs=-1).fit(X_pca)
    distances, indices = nn_res.kneighbors(X_pca)
    G = nx.Graph()
    for i in range(indices.shape[0]):
        for j in indices[i][1:]:
            G.add_edge(i, j)
    partition = community_louvain.best_partition(G, random_state=seed)
    return partition

def write_list_to_file(lst, filename):
    with open(filename, 'w') as file:
        for sublist in lst:
            line = '\t'.join(str(element) for element in sublist)
            file.write(line + '\n')

"""
scCAD: Cluster decomposition-based Anomaly Detection for Rare Cell Identification in Single Cell Expression Data.
"""
def scCAD(data,
          dataName=None,
          cellNames=None,
          geneNames=None,
          normalization=True,
          seed=2023,
          merge_h=50,
          overlap_h=0.7,
          rare_h=0.01,
          save_full=True,
          save_path='./'):
    """
    Parameters
    ----------
    data : `pandas.DataFrame` or `2-D numpy.array`, optional
        The expression data matrix with genes in columns and cells in rows.

    dataName : string
        Name of scRNA-seq dataset. (default: None)

    cellNames : list -> string
        The length must be the same as the number of rows in the data matrix.
        Names of all cells. (default: None)

    geneNames : list -> string
        The length must be the same as the number of columns in the data matrix.
        Names of all genes. (default: None)

    normalization : boolean
        Whether the data needs to be normalized. (default: True)

    seed : integer
        Random seed. (default: 2023)

    merge_h : float
        Threshold to use when doing cluster merge. (default: 50)

    overlap_h : float
        Overlap threshold to identify rare clusters. (default: 0.7)

    rare_h : float
        Rare threshold to use when doing cluster decomposition. (default: 0.01)

    save_full : boolean
        Whether the full result needs to be saved. (default: True)

    save_path : string
        Path to save results.

    Returns
    ----------
    result: `list`
        rare cell list

    remain_degs_list: `list`
        rare cell de genes list

    rename_comb_subclusters: 'np.array'
        sub-clusters assignment for all cells, length: n cells

    overlap: 'np.array'
        independence scores assignment for all sub-clusters, length: n sub-clusters

    """

    warnings.filterwarnings("ignore")
    if dataName is None:
        dataName = ""

    # 0. Load data
    adata = sc.AnnData(data)
    if geneNames is not None:
        geneNames = np.array(geneNames)
        adata.var_names = geneNames
        if len(geneNames) != len(np.unique(geneNames)):
            print(">>> The gene name is repeated and changed to a unique identifier!")
            adata.var_names = [str(i) for i in range(adata.X.shape[1])]
    else:
        adata.var_names = [str(i) for i in range(adata.X.shape[1])]
    if cellNames is not None:
        cellNames = np.array(cellNames)
        adata.obs_names = cellNames
    else:
        adata.obs_names = [str(i) for i in range(adata.X.shape[0])]

    # 1. Data preprocessing
    if normalization:
        print(">>> Data preprocessing in progress...")
        adata = normalize(adata)

    X_norm_withallgenes = adata.X.copy()
    n_cells = X_norm_withallgenes.shape[0]
    n_genes = X_norm_withallgenes.shape[1]
    if normalization:
        print(">>> After preprocessing, cells: %d; genes: %d;" % (n_cells, n_genes))
    else:
        print(">>> Cells: %d; genes: %d;" % (n_cells, n_genes))

    start = time.time()

    # 2 Feature selection
    print(">>> Feature selection in progress...")
    # 2.1 HVGs
    print(">>> Identification of HVGs is currently in progress...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    hvgs = list(adata[:, adata.var.highly_variable].var_names)

    # 2.2 RFGs
    print(">>> Identification of RFGs is currently in progress...")
    Init_clusters = fast_clustering(data=X_norm_withallgenes, seed=seed)
    Init_clusters = [Init_clusters[node] for node in range(n_cells)]
    rf = RandomForestClassifier(n_estimators=1000, n_jobs=-1, random_state=seed, bootstrap=False)
    rf.fit(X_norm_withallgenes, Init_clusters)
    gene_imp = rf.feature_importances_.copy()
    rfgs = list(adata.var_names[np.argsort(-gene_imp)[:2000]])
    if save_full:
        np.savetxt(save_path + dataName + '_scCAD_Init_clusters.txt', Init_clusters, fmt='%d')
    sg = list(set(hvgs).union(set(rfgs)))

    adata = adata[:, sg]
    X_norm = adata.X.copy()
    print(">>> After feature selection, genes: %d;" % (X_norm.shape[1]))
    if save_full:
        np.savetxt(save_path + dataName + '_scCAD_select_genes.txt', sg, fmt='%s')

    # 3 Clusters decomposition
    print(">>> Clusters decomposition in progress...")
    pseudo_init_subclusters = fast_clustering(data=X_norm, seed=seed)
    pseudo_init_subclusters = np.array([pseudo_init_subclusters[node] for node in range(n_cells)])
    if save_full:
        np.savetxt(save_path + dataName + '_scCAD_Init_clusters_based_selected_genes.txt', pseudo_init_subclusters, fmt='%d')

    pseudo_subclusters = pseudo_init_subclusters.copy()
    h1 = max(int(rare_h * n_cells), 30)
    count = 1
    while 1:
        print(">>> iter %d, running..." % count)
        count = count + 1
        dict = Counter(pseudo_subclusters)
        depths = {}
        dpt = 1
        while dict.most_common(1)[0][1] >= h1:
            depths.update(
                (key, dpt) for key in list(set([i for i, count in dict.items() if count < h1]) - set(depths.keys())))
            dpt = dpt + 1
            c_max = max(pseudo_subclusters) + 1
            c_list = list(set([i for i, count in dict.items() if count >= h1]) - set(depths.keys()))
            if len(c_list) == 0:
                break
            for clustid in c_list:
                idx = np.where(pseudo_subclusters == clustid)[0]
                temp_X = X_norm[idx, :].copy()
                temp_clusters = fast_clustering(data=temp_X, seed=seed)
                temp_clusters = [temp_clusters[node] + c_max for node in range(temp_X.shape[0])]
                if len(np.unique(temp_clusters)) != 1:
                    pseudo_subclusters[idx] = temp_clusters
                    c_max = max(pseudo_subclusters) + 1
                else:
                    depths[clustid] = dpt - 1
            dict = Counter(pseudo_subclusters)

        subc_dict = {}
        p = 0
        rename_pseudo_subclusters = np.zeros(n_cells, dtype=int)
        for i in range(n_cells):
            if pseudo_subclusters[i] in subc_dict.keys():
                rename_pseudo_subclusters[i] = subc_dict[pseudo_subclusters[i]]
            else:
                subc_dict[pseudo_subclusters[i]] = p
                rename_pseudo_subclusters[i] = p
                p = p + 1
        n_subclusters = len(np.unique(rename_pseudo_subclusters))
        print(">>> After clusters decomposition, we got %d balanced sub-clusters." % (n_subclusters))

        # 4 Clusters merge
        print(">>> Clusters merge in progress...")
        X_centers = np.zeros((len(np.unique(rename_pseudo_subclusters)), X_norm.shape[1]))
        for i in np.unique(rename_pseudo_subclusters):
            id = np.where(rename_pseudo_subclusters == i)[0]
            X_centers[i, :] = np.mean(X_norm[id, :], axis=0)

        distances = pdist(X_centers)
        distance_matrix = squareform(distances)
        np.fill_diagonal(distance_matrix, np.max(distance_matrix) + 1)
        comb_id = []
        for i in range(X_centers.shape[0]):
            i_nearest_neighbor_index = np.argsort(distance_matrix[i, :])
            for j in i_nearest_neighbor_index:
                dist = distance_matrix[i, j]
                if dist <= np.percentile(np.min(distance_matrix, axis=0), merge_h) and sorted([i, j]) not in comb_id:
                    comb_id.append(sorted([i, j]))

        merged = []
        for sublist in comb_id:
            merged_with = None
            for m in merged:
                if any(x in m for x in sublist):
                    merged_with = m
                    break
            if merged_with:
                merged.remove(merged_with)
                merged.append(list(set(merged_with + sublist)))
            else:
                merged.append(sublist)

        comb_subclusters = rename_pseudo_subclusters.copy()
        for i in tqdm(range(len(merged))):
            for j in range(n_cells):
                if comb_subclusters[j] in merged[i]:
                    comb_subclusters[j] = min(merged[i])

        subc_dict = {}
        p = 0
        rename_comb_subclusters = np.zeros(n_cells, dtype=int)
        for i in range(n_cells):
            if comb_subclusters[i] in subc_dict.keys():
                rename_comb_subclusters[i] = subc_dict[comb_subclusters[i]]
            else:
                subc_dict[comb_subclusters[i]] = p
                rename_comb_subclusters[i] = p
                p = p + 1
        n_subclusters = len(np.unique(comb_subclusters))
        print(">>> After clusters merge, we got %d sub-clusters." % (n_subclusters))

        # 5 Cluster anomaly score calculation
        print(">>> Cluster anomaly score calculation in progress...")
        IFmodel = IsolationForest(n_estimators=100, random_state=seed, n_jobs=-1)
        overlap = []
        degs_list = []
        for i in tqdm(range(n_subclusters)):
            id = np.where(rename_comb_subclusters == i)[0]
            if len(id) > h1 or len(id) < 10:
                overlap.append(0)
                degs_list.append([])
            else:
                id_ = np.where(rename_comb_subclusters != i)[0]
                tmp_X = X_norm[id, :]
                zero_cols = np.where(np.all(tmp_X == 0, axis=0))[0]
                re_cols = list(set(np.arange(X_norm.shape[1])) - set(zero_cols))
                tmp_X = X_norm[:, re_cols]
                diff = np.abs(np.median(tmp_X[id, :], axis=0) - np.median(tmp_X[id_, :], axis=0))
                var_names = np.array(sg)[re_cols]
                n_top = min(20, len(np.where(diff > 0)[0]))
                degs_ = list(var_names[np.argsort(-diff)[:n_top]])
                degs_list.append(degs_)
                IFmodel.fit(adata[:, degs_].X)
                s = IFmodel.score_samples(adata[:, degs_].X)
                overlap.append(len(set(np.argsort(s)[:len(id)]) & set(id)) / len(id))

        remain_clusters = []
        overlap = np.array(overlap)
        remain_degs_list = []
        for i in range(n_subclusters):
            if overlap[i] >= overlap_h:
                remain_clusters.append(i)
                remain_degs_list.append(degs_list[i])

        if len(remain_clusters) != 0:
            break
        elif rare_h > 0.05:
            print(">>> Rare type (<5%) not found!")
            break
        else:
            if 30 >= int(rare_h * n_cells):
                rare_h = 30/n_cells + 0.01
            else:
                rare_h = rare_h + 0.01
            h1 = int(rare_h * n_cells)
            pseudo_subclusters = pseudo_init_subclusters.copy()

    end = time.time()
    runningtime = end - start
    print(">>> time used:", runningtime)

    # 6 Output results
    result = []
    for i in remain_clusters:
        id = np.where(rename_comb_subclusters == i)[0]
        if cellNames is None:
            result.append(list(id))
        else:
            result.append(list(cellNames[id]))

    if save_full:
        np.savetxt(save_path + dataName + '_scCAD_balanced_sub-clusters.txt', rename_pseudo_subclusters, fmt='%d')
        np.savetxt(save_path + dataName + '_scCAD_comb_sub-clusters.txt', rename_comb_subclusters, fmt='%d')
        np.savetxt(save_path + dataName + '_scCAD_sub-clusters_anomaly_score.txt', overlap, fmt='%f')

    write_list_to_file(result, save_path + dataName + '_scCAD_rare_cells_result.txt')
    write_list_to_file(remain_degs_list, save_path + dataName + '_scCAD_degs_list.txt')

    return result, overlap, rename_comb_subclusters, remain_degs_list

