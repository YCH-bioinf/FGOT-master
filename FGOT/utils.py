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

def graph_alpha(spatial_locs, n_neighbors=20):
    """
    Construct a geometry-aware spatial proximity graph of the spatial spots of cells by using alpha complex.
    """
    A_knn = kneighbors_graph(spatial_locs, n_neighbors=n_neighbors, mode='distance')
    estimated_graph_cut = A_knn.sum() / float(A_knn.count_nonzero())
    spatial_locs_list = spatial_locs.tolist()
    n_node = len(spatial_locs_list)
    alpha_complex = gudhi.AlphaComplex(points=spatial_locs_list)
    simplex_tree = alpha_complex.create_simplex_tree(max_alpha_square=estimated_graph_cut ** 2)
    skeleton = simplex_tree.get_skeleton(1)
    initial_graph = nx.Graph()
    initial_graph.add_nodes_from([i for i in range(n_node)])
    for s in skeleton:
        if len(s[0]) == 2:
            initial_graph.add_edge(s[0][0], s[0][1])
    extended_graph = nx.Graph()
    extended_graph.add_nodes_from(initial_graph)
    extended_graph.add_edges_from(initial_graph.edges)
    return nx.to_scipy_sparse_array(extended_graph, format='csr')

def graph_high_gene_expression(extended_graph, gene_expression,ratio=0.75):
    """
    gene_expression: numpy.ndarray
    extended_graph: csr scipy_sparse_array
    """
    if gene_expression.shape[0] != extended_graph.shape[0]:
        print("the length of cell is not equal of the length of spatial locations!")
    for i in range(extended_graph.shape[0]):
        _, row_cols = extended_graph[[i]].nonzero()
        similarities = cosine_similarity([gene_expression[i,]], gene_expression[row_cols,])[0]
        sorted_indices = np.argsort(similarities)[::-1]
        top_75_percent = int(ratio * len(similarities))
        top_25_percent_indices = sorted_indices[top_75_percent:]
        for j in top_25_percent_indices:
            col = row_cols[j]
            extended_graph[i,col] = 0
    return extended_graph

def graph_imagefeature(image_features, percentile=75):
    """
    Construct a graph using image features.
    """
    graph = nx.Graph()
    distances = np.linalg.norm(image_features[:, np.newaxis] - image_features, axis=2)
    threshold = np.percentile(distances, percentile)
    indices = np.where(distances <= threshold)
    edges = list(zip(indices[0], indices[1]))
    graph.add_edges_from(edges)
    return nx.to_scipy_sparse_array(graph, format='csr')

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

def select_confident_cells(adata, label_name, q_cutoff = 0.75,focus_cell_type=None):
    adata_sc = adata.copy()
    sc.pp.scale(adata_sc)
    # calculate leiden clustering
    sc.pp.pca(adata_sc)
    sc.pp.neighbors(adata_sc)
    snn_mat = np.asmatrix(adata_sc.obsp['connectivities'].A)
    
    label_type = adata_sc.obs[label_name].tolist()
    if focus_cell_type is None:
        focus_cell_type = list(set(label_type))
    
    min_cluster = 20
    index_type_l = []
    for i in range(len(focus_cell_type)):
        indexes = [index for index in range(len(label_type)) if label_type[index] == focus_cell_type[i]]
        if len(indexes) > min_cluster:
            snn_i = snn_mat[np.ix_(indexes, indexes)]
            col_si = snn_i.mean(1)
            thresh_i = np.quantile(np.array(col_si)[:,0], q_cutoff)
            index_i = np.where(col_si > thresh_i)[0]
            if len(index_i) > min_cluster:
                index_type_l.extend([indexes[index_i[j]] for j in range(len(index_i))])
            else:
                index_type_l.extend([indexes[j] for j in range(len(indexes))])
        else:
            index_type_l.extend([indexes[j] for j in range(len(indexes))])
    anchor_indnums=list(set(index_type_l))
    adata = adata[anchor_indnums,:]
    anchor_names = list(adata.obs.index)
    return anchor_names,anchor_indnums

def compute_D_kl(Snn1, Snn2, k, l)->np.ndarray:
    row_k = Snn1[k, :]
    row_l = Snn2[l, :]
    row_k = row_k[:, np.newaxis]
    row_l = row_l[np.newaxis, :]
    difference_square = (row_k - row_l) ** 2 # Numpy broadcasting
    D_kl = 0.5 * difference_square
    
    return D_kl
