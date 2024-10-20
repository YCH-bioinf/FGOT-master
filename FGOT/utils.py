import numpy as np
from sklearn.neighbors import NearestNeighbors
import hnswlib
from sklearn.neighbors import kneighbors_graph
import gudhi
import networkx as nx
from scipy import sparse
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import pdist, squareform
import scanpy as sc
from anndata import AnnData
import torch


def jaccard_index(A,B):
    set_A, set_B = set(A),set(B)
    intersection_size = len(set_A.intersection(set_B))
    union_size = len(set_A.union(set_B))
    jaccard_index = intersection_size / union_size
    return jaccard_index

def nn_approx(ds1, ds2, names1, names2, knn=50):
    dim = ds2.shape[1]
    num_elements = ds2.shape[0]
    p = hnswlib.Index(space='l2', dim=dim)
    p.init_index(max_elements=num_elements, ef_construction=100, M = 16)
    p.set_ef(10)
    p.add_items(ds2)
    ind,  distances = p.knn_query(ds1, k=knn)
    match = set()
    for a, b in zip(range(ds1.shape[0]), ind):
        for b_i in b:
            match.add((names1[a], names2[b_i]))
    return match

def nn(ds1, ds2, names1, names2, knn=50):
    # Find nearest neighbors of first dataset.
    nn_ = NearestNeighbors(n_neighbors= knn)
    nn_.fit(ds2)
    ind = nn_.kneighbors(ds1, return_distance=False)

    match = set()
    for a, b in zip(range(ds1.shape[0]), ind):
        for b_i in b:
            match.add((names1[a], names2[b_i]))
    return match

def mnn(ds1, ds2, names1, names2, knn = 20, save_on_disk = True, approx = True):
    if approx: 
        # Find nearest neighbors in first direction.
        # output KNN point for each point in ds1.  match1 is a set(): (points in names1, points in names2), the size of the set is ds1.shape[0]*knn
        match1 = nn_approx(ds1, ds2, names1, names2, knn=knn)
        # Find nearest neighbors in second direction.
        match2 = nn_approx(ds2, ds1, names2, names1, knn=knn)
    else:
        match1 = nn(ds1, ds2, names1, names2, knn=knn)
        match2 = nn(ds2, ds1, names2, names1, knn=knn)
    # Compute mutual nearest neighbors.
    mutual = match1 & set([ (b, a) for a, b in match2 ])
    return mutual

def get_distance_matrix(spatial):
    spot_coordinates = spatial
    distances = pdist(spot_coordinates)
    distance_matrix = squareform(distances)
    return distance_matrix

def coo_submatrix_pull(matr, rows, cols):
    """
    Pulls out an arbitrary i.e. non-contiguous submatrix out of
    a sparse.coo_matrix. 
    """
    if type(matr) != sparse.coo_matrix:
        raise TypeError('Matrix must be sparse COOrdinate format')
    
    gr = -1 * np.ones(matr.shape[0])
    gc = -1 * np.ones(matr.shape[1])
    
    lr = len(rows)
    lc = len(cols)
    
    ar = np.arange(0, lr)
    ac = np.arange(0, lc)
    gr[rows[ar]] = ar
    gc[cols[ac]] = ac
    mrow = matr.row
    mcol = matr.col
    newelem = (gr[mrow] > -1) & (gc[mcol] > -1)
    newrows = mrow[newelem]
    newcols = mcol[newelem]
    return sparse.coo_matrix((matr.data[newelem], np.array([gr[newrows],
        gc[newcols]])),(lr, lc))

def fgot_flatten_matrix_sparse(M:np.ndarray, S:np.ndarray, D:np.ndarray, A:np.ndarray, ns_s:int, ns_d:int, n_pos_s:int, n_pos_d:int):
    a = S.flatten('F')
    b = D.flatten('F')
    
    C_data, C_row, C_col = [], [], []
    M_row, M_col = np.indices(M.shape)
    M_row, M_col = M_row.reshape(-1), M_col.reshape(-1)
    M_max_sp = sparse.coo_matrix((M[M_row,M_col], (M_row,M_col)), shape=M.shape)
    del M_row, M_col
    
    cost_scales = []
    for i in range(ns_s):
        for j in range(ns_d):
            if not np.isinf(A[i,j]):
                tmp_nzind_s = np.where(S[:,i] > 0)[0]
                tmp_nzind_d = np.where(D[:,j] > 0)[0]
                tmp_M_max_sp = coo_submatrix_pull(M_max_sp, tmp_nzind_s, tmp_nzind_d)
                tmp_row = tmp_nzind_s[tmp_M_max_sp.row]
                tmp_col = tmp_nzind_d[tmp_M_max_sp.col]
                cost_scales.append(np.max(M_max_sp.data) * A[i, j])
                C_data.append( tmp_M_max_sp.data*A[i,j] )
                C_row.append( tmp_row+i*n_pos_s )
                C_col.append( tmp_col+j*n_pos_d )
    del S, D,tmp_nzind_s,tmp_nzind_d,tmp_M_max_sp,tmp_row,tmp_col
    cost_scale = np.max(cost_scales)
    del cost_scales
    C_data = np.concatenate(C_data, axis=0)
    C_row = np.concatenate(C_row, axis=0)
    C_col = np.concatenate(C_col, axis=0)
    C = sparse.coo_matrix((C_data/cost_scale, (C_row, C_col)), shape=(len(a),len(b)))
    
    return C
    
def get_snn_mat(X:np.ndarray, scale = 1, n_neighbors = 50):
    n = X.shape[0]
    andata = AnnData(X)
    if scale == 1:
        sc.pp.scale(andata)
    sc.tl.pca(andata, n_comps = 50)
    sc.pp.neighbors(andata, n_neighbors= n_neighbors, n_pcs=15)
    snn = andata.obsp['connectivities'].todense()
    snn = np.array(snn)
    snn[range(n), range(n)] = 1
    return snn

def select_central_cells(S, label_name, q_cutoff = 0.75, focus_cell_type=None):
    if focus_cell_type is None:
        focus_cell_type = list(set(label_name))
    index_type_l = []
    S = np.array(S)
    for i in range(len(focus_cell_type)):
        indexes = [index for index in range(len(label_name)) if label_name[index] == focus_cell_type[i]]
        wnn_i = S[np.ix_(indexes, indexes)]
        col_si = wnn_i.mean(axis=0)
        thresh_i = np.quantile(np.array(col_si), q_cutoff)
        index_i = np.where(col_si > thresh_i)[0]
        index_type_l.extend([indexes[j] for j in index_i])
    anchor_indnums = list(set(index_type_l))
    return anchor_indnums


def select_anchor_cells_unpaired(S, D, mnn_df, label_name1, label_name2, q_cutoff = 0.75, focus_cell_type=None, scale = 0):
    snn1 = get_snn_mat(S, scale=scale)
    snn2 = get_snn_mat(D, scale=scale)
    index_central_cellsl = select_central_cells(snn1, label_name1, q_cutoff = q_cutoff, focus_cell_type=focus_cell_type)
    index_central_cells2 = select_central_cells(snn2, label_name2, q_cutoff = q_cutoff, focus_cell_type=focus_cell_type)
    mnn_pairs_indices = np.argwhere(mnn_df.values == 1)
    anchor_indnumsS = []
    anchor_indnumsD = []
    for row, col in mnn_pairs_indices:
        if row in index_central_cellsl or col in index_central_cells2:
            anchor_indnumsS.append(row)
            anchor_indnumsD.append(col)
    anchor_indnumsS = list(set(anchor_indnumsS))
    anchor_indnumsD = list(set(anchor_indnumsD))
    return anchor_indnumsS, anchor_indnumsD


def compute_cost_anchor(Snn1_batch, Snn2_batch, P_anchor, device='cuda'):
    Snn1_batch = torch.tensor(Snn1_batch, device=device, dtype=torch.float32)
    Snn2_batch = torch.tensor(Snn2_batch, device=device, dtype=torch.float32)
    P_anchor = torch.tensor(P_anchor, device=device, dtype=torch.float32)

    m = Snn1_batch.shape[0]
    n = Snn2_batch.shape[0]
    vector_m = torch.ones((m, 1), device=device, dtype=torch.float32)
    vector_n = torch.ones((n, 1), device=device, dtype=torch.float32)
    p = torch.sum(P_anchor, axis=1, keepdim=True)
    q = torch.sum(P_anchor, axis=0, keepdim=True).T
    C = (Snn1_batch ** 2) @ p @ vector_n.T + vector_m @ q.T @ (Snn2_batch ** 2).T - Snn1_batch @ P_anchor @ Snn2_batch.T

    return C.cpu().numpy() if device != 'cpu' else C.numpy()