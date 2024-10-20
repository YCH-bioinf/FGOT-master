import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, normalize
import networkx as nx
import torch
from tqdm import tqdm

from .utils import jaccard_index,mnn,get_snn_mat

def calculate_group_affinity(markers1,markers2):
    batch1_label2marker = {}
    for gene in list(markers1['gene']):
        group = markers1.loc[gene].cluster
        if not group in batch1_label2marker:
            batch1_label2marker[group] = []
        batch1_label2marker[group].append(gene)

    batch2_label2marker = {}
    for gene in list(markers2['gene']):
        group = markers2.loc[gene].cluster
        if not group in batch2_label2marker:
            batch2_label2marker[group] = []
        batch2_label2marker[group].append(gene)
        
    group_affinity = pd.DataFrame(0, index=batch1_label2marker.keys(), columns=batch2_label2marker.keys())
    for group_i in group_affinity.index:
        for group_j in group_affinity.columns:
            group_affinity.loc[group_i,group_j] = jaccard_index(batch1_label2marker[group_i],batch2_label2marker[group_j])
    return group_affinity

def calculate_cell_similarity_byMNN(X1, X2, cell_names1, cell_names2, scale=1, knn=50):
    # scale and then calculate MNN pairs
    if scale == 1:
        scaler = StandardScaler()
        X1, X2 = scaler.fit_transform(X1), scaler.fit_transform(X2)
    cell_names2 = [str+"-1" for str in cell_names2]
    
    match = mnn(X1, X2, cell_names1, cell_names2, knn=knn, save_on_disk = False, approx = False)
    print("the number of MNN pairs is ",len(match))
    G = nx.Graph()
    nodes = np.hstack((cell_names1, cell_names2))
    G.add_nodes_from(nodes)
    G.add_edges_from(match)
    adj = nx.adjacency_matrix(G)
    n1 = len(cell_names1)
    mnn_mtx = adj.todense()[:n1,n1:]
    mnn_mtx = pd.DataFrame(mnn_mtx, index = cell_names1, columns=cell_names2)
    return mnn_mtx

# check the ratio of accuracy
def check_mnn_accuracy(mnn_adj, label1, label2):
    mnn_adj = np.array(mnn_adj)
    con_same = np.zeros((mnn_adj.shape[0],2))
    for p in range(mnn_adj.shape[0]):
        c_p = mnn_adj[p,:].nonzero()[0]
        if len(c_p) > 0:
            label_cp = label2[c_p]
            con_same[p,0] = sum(label_cp == label1[p])/len(c_p)
            con_same[p,1] = 1
    r = sum(con_same[con_same[:,1].nonzero(),0][0])/len(con_same[:,1].nonzero()[0])
    return r

def smooth_cell_similarity_byLaplacian1(cell_Cor, X1, X2, n_neighbors = 50, eps = 1e-12):
    snn1 = get_snn_mat(X1,n_neighbors)
    snn2 = get_snn_mat(X2,n_neighbors)
    A1 = torch.tensor(snn1, dtype = torch.float32)
    A2 = torch.tensor(snn2, dtype = torch.float32)
    D1 = torch.diag(A1.sum(1))
    D2 = torch.diag(A2.sum(1))
    deg1 = A1.sum(dim = 1)
    deg2 = A2.sum(dim = 1)
    deg1 += eps
    deg2 += eps
    D1_sqrt_inv = 1/deg1.sqrt()
    D1_sqrt_inv = torch.diag(D1_sqrt_inv)
    D2_sqrt_inv = 1/deg2.sqrt()
    D2_sqrt_inv = torch.diag(D2_sqrt_inv)
    L1 = D1_sqrt_inv@A1@D1_sqrt_inv
    L2 = D2_sqrt_inv@A2@D2_sqrt_inv
    M = torch.tensor(np.array(cell_Cor), dtype = torch.float32)
    S = L1@M@L2
    S_new = pd.DataFrame(S, index = cell_Cor.index, columns=cell_Cor.columns)
    return S_new

def smooth_cell_similarity_byLaplacian2(cell_Cor, snn1, snn2, eps = 1e-12):
    snn1 = np.array(snn1)
    snn2 = np.array(snn2)
    A1 = torch.tensor(snn1, dtype = torch.float32)
    A2 = torch.tensor(snn2, dtype = torch.float32)
    D1 = torch.diag(A1.sum(1))
    D2 = torch.diag(A2.sum(1))
    deg1 = A1.sum(dim = 1)
    deg2 = A2.sum(dim = 1)
    deg1 += eps
    deg2 += eps
    D1_sqrt_inv = 1/deg1.sqrt()
    D1_sqrt_inv = torch.diag(D1_sqrt_inv)
    D2_sqrt_inv = 1/deg2.sqrt()
    D2_sqrt_inv = torch.diag(D2_sqrt_inv)
    L1 = D1_sqrt_inv@A1@D1_sqrt_inv
    L2 = D2_sqrt_inv@A2@D2_sqrt_inv
    M = torch.tensor(np.array(cell_Cor), dtype = torch.float32)
    S = L1@M@L2
    S_new = pd.DataFrame(S, index = cell_Cor.index, columns=cell_Cor.columns)
    return S_new

def prior_feature_graph(promoters, peak_names, gene_names, scope = 250000):
    filtered_promoters = []
    columns = ['id', 'chr', 'starts', 'ends', 'forward', 'backward', 'gene']
    for promoter in tqdm(promoters.itertuples()):
        chr = promoter.chr
        starts = promoter.starts
        ends = int(promoter.ends)
        genes = promoter.genes
        promoter_id = chr + '-' + str(starts) + '-' + str(ends) + '-' + genes
        forward = max(0, starts - scope)
        backward = starts + scope
        filtered_promoters.append([promoter_id, chr, starts, ends, forward, backward, genes])
    filtered_promoters = pd.DataFrame(filtered_promoters, columns=columns)
    filtered_promoters = filtered_promoters.drop_duplicates('id')
    gene_peaks = {}
    for promoter in tqdm(filtered_promoters.itertuples()):
        if not promoter.gene in gene_names:
            continue
        id = promoter.id
        chr = promoter.chr
        starts = promoter.starts
        ends = promoter.ends
        forward = promoter.forward
        backward = promoter.backward
        gene = promoter.gene
        if not gene in gene_peaks:
            gene_peaks[gene] = set()
        for peak in peak_names:
            peak_chr, peak_start, peak_end = peak.split('-')
            if peak_chr == chr and int(peak_start) >= forward and int(peak_end) <= backward:
                gene_peaks[gene].add(peak)
    for key in gene_peaks.keys():
        gene_peaks[key] = list(gene_peaks[key])
    gene_nodes = []
    peak_nodes = []
    edges = []
    for key in tqdm(gene_peaks.keys()):
        if len(gene_peaks[key]) == 0:
            continue
        if key not in gene_nodes:
            gene_nodes.append(key)
        peaks = gene_peaks[key]
        for peak in peaks:
            if peak in peak_names:
                if peak not in peak_nodes:
                    peak_nodes.append(peak)
                edge = [key, peak]
                edges.append(edge)
    print("The number of gene nodes, peak nodes, and edges in the prior feature graph is:",len(gene_nodes),len(peak_nodes),len(edges))
    hvg_indexs = {gene_names[i]:i for i in range(len(gene_names))}
    peak_indexs = {peak_names[i]:i for i in range(len(peak_names))}
    adj = np.zeros((len(peak_names), len(gene_names)))
    for edge in tqdm(edges):
        hvg_index = hvg_indexs[edge[0]]
        peak_index = peak_indexs[edge[1]]
        adj[peak_index, hvg_index] = 1
    adj = pd.DataFrame(adj, index=peak_names, columns=gene_names)
    adj[adj == 0] = np.inf
    return adj
    