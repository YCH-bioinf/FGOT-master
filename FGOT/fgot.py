import numpy as np
import pandas as pd
from scipy import sparse
import random
from tqdm import tqdm
import torch
from sklearn.neighbors import NearestNeighbors
import scanpy as sc
from anndata import AnnData
from .fgot_solver import fgot_flatten_solve
from .utils import select_central_cells, select_anchor_cells_unpaired, compute_cost_anchor, get_snn_mat
from collections import defaultdict


def fgot_tol(P_tensor):
    shape = P_tensor.get(list(P_tensor.keys())[0]).shape
    P_2dim = np.zeros(shape)
    for key in tqdm(P_tensor.keys()):
        sparse_matrix = P_tensor[key].todense()
        P_2dim += sparse_matrix
    return P_2dim

def fgot_sparse_tensor(S:pd.DataFrame, D:pd.DataFrame, A : pd.DataFrame, M: np.ndarray, \
                    cell1_cluster=None, cell2_cluster=None, minibatch=0, batchsize=400, pair = False,\
                    eps_p=1e-1, rho_mu=1e1, rho_nu=1e1, nitermax=1e3, stopthr=1e-8, device='cpu',\
                    fastMinibatch=True, P_anchor:pd.DataFrame=None, lam:float=0.5,adjust=True)->dict:
    """ adjust
    Set batches for mini-batch computing.
    Set various parameters for feature-guided optimal transport problem.    
    Parameters
    ----------
    S : (n_pos_d,ns_d) pandas.dataframe
        Source distributions over `n_pos_s` positions of `ns_s` source species.
    D : (n_pos_d,ns_d) pandas.dataframe
        Destination distributions over `n_pos_d` positions of `ns_d` destination species.
    A : (ns_s,ns_d) pandas.dataframe
        The cost coefficients for source-destination species pairs. An infinity value indicates that the two species cannot be coupled.
        feature matrix
    M : (n_pos_s,n_pos_d) pandas.dataframe
        The distance (cost) matrix among the positions.
        cost
    
    minbatch: int
        0 or 1
    Returns
    -------
    A dictionary of scipy.sparse.coo_matrix
        The transport plan in coo sparse format for source species i and destinaton species j can be retrieved with the key (i,j).
    """ 
    if minibatch == 0:
        print("minibatch = 0 and device = ",device)
        P_cot = fgot_flatten_solve(S, D, A, M,\
            eps_p=eps_p, rho_mu=rho_mu, rho_nu=rho_nu,\
            nitermax=nitermax, stopthr=stopthr, device=device)
    if minibatch == 1:
        if cell1_cluster is None or cell2_cluster is None:
            print("Please input both cell1_cluster and cell2_cluster!")
            return
        print("minibatch = 1 and device = ",device)
        n1 = len(S.index)
        n2 = len(D.index)
        samples1 = sample_from_cluster(cell1_cluster,batchsize)
        if pair:
            if list(S.index) != list(D.index):
                print("Please ensure that the data is pair!")
                return
            samples2 = samples1
            total_batch = len(samples1)
        else:
            samples2 = sample_from_cluster(cell2_cluster,batchsize)
            total_batch = len(samples1) * len(samples2)
        print("the total number of computing batch is ", total_batch)
        
        if fastMinibatch==False: # precise minibatch
            if(P_anchor is None):
                print('Please input P_anchor!')
            anchors1 = P_anchor.index
            anchors2 = P_anchor.columns
            anchors1_id = [S.index.get_loc(x) for x in anchors1]
            anchors2_id = [D.index.get_loc(x) for x in anchors2]
            P_anchor = P_anchor.values
            Snn1 = get_snn_mat(S, scale=0)
            Snn2 = get_snn_mat(D, scale=0)
        P_cot = {}
        pbar = tqdm(total=total_batch)
        # tmp_cnt = 0 # for test output
        
        for cells1 in samples1:
            if pair:
                cells2 = cells1
                if len(cells1) == 0 or len(cells2) == 0:
                    print("the set of sample is not right!")
                    return
                cells1_id = np.array(cells1)
                cells2_id = np.array(cells2)
                S_batch = S.iloc[cells1_id]
                D_batch = D.iloc[cells2_id]
                M_batch = M.iloc[cells1_id,cells2_id]
                
                if adjust:
                    r_sampleRate = len(cells1_id) / n1
                    c_sampleRate = len(cells2_id) / n2
                    rho_mu = rho_mu + eps_p * np.log(1/c_sampleRate)
                    rho_nu = rho_nu + eps_p * np.log(1/r_sampleRate)
                else:
                    r_sampleRate = 1
                    c_sampleRate = 1
                
                if fastMinibatch==False:
                    Snn1_batch = Snn1[cells1_id,:][:,anchors1_id]
                    Snn2_batch = Snn2[cells2_id,:][:,anchors2_id]
                    P_anchor = np.where(P_anchor < 1e-4, 0, P_anchor)
                    M_anchor = compute_cost_anchor(Snn1_batch, Snn2_batch, P_anchor, device=device)
                    torch.cuda.empty_cache()
                    lam = lam * (np.max(M_batch)/np.max(M_anchor))
                    M_batch = M_batch + lam * M_anchor
                    
                P_tmp = fgot_flatten_solve(S_batch, D_batch, A, M_batch,\
                    eps_p=eps_p, rho_mu=rho_mu, rho_nu=rho_nu,\
                    nitermax=nitermax, stopthr=stopthr,device=device,\
                    c_sampleRate=c_sampleRate,r_sampleRate=r_sampleRate)
                torch.cuda.empty_cache()
                
                for key in P_tmp.keys():
                    sparse_mtx = P_tmp.get(key)
                    P_tmp[key] = sparse.coo_matrix((sparse_mtx.data, (cells1_id[sparse_mtx.row], cells2_id[sparse_mtx.col])), shape=(n1,n2))
                    if not (key in P_cot.keys()):
                        P_cot[key]= P_tmp[key]
                    else:
                        P_cot[key].data = np.concatenate((P_cot[key].data,P_tmp[key].data))
                        P_cot[key].row  = np.concatenate((P_cot[key].row,P_tmp[key].row))
                        P_cot[key].col  = np.concatenate((P_cot[key].col,P_tmp[key].col))
                
                pbar.update()
            else:
                for cells2 in samples2:   
                    if len(cells1) == 0 or len(cells2) == 0:
                        print("the set of sample is not right!")
                        return
                    cells1_id = np.array(cells1)
                    cells2_id = np.array(cells2)
                    S_batch = S.iloc[cells1_id]
                    D_batch = D.iloc[cells2_id]
                    M_batch = M.iloc[cells1_id,cells2_id]
                    
                    if adjust:
                        r_sampleRate = len(cells1_id) / n1
                        c_sampleRate = len(cells2_id) / n2
                        rho_mu = rho_mu + eps_p * np.log(1/c_sampleRate)
                        rho_nu = rho_nu + eps_p * np.log(1/r_sampleRate)
                    else:
                        r_sampleRate = 1
                        c_sampleRate = 1

                    if fastMinibatch==False:
                        Snn1_batch = Snn1[cells1_id,:][:,anchors1_id]
                        Snn2_batch = Snn2[cells2_id,:][:,anchors2_id]
                        P_anchor = np.where(P_anchor < 1e-4, 0, P_anchor)
                        M_anchor = compute_cost_anchor(Snn1_batch, Snn2_batch, P_anchor, device=device)
                        torch.cuda.empty_cache()
                        lam = lam * (np.max(M_batch)/np.max(M_anchor))
                        M_batch = M_batch + lam * M_anchor
                    
                    P_tmp = fgot_flatten_solve(S_batch, D_batch, A, M_batch,\
                        eps_p=eps_p, rho_mu=rho_mu, rho_nu=rho_nu,\
                        nitermax=nitermax, stopthr=stopthr,device=device,\
                        c_sampleRate=c_sampleRate,r_sampleRate=r_sampleRate)
                    torch.cuda.empty_cache()
                    
                    for key in P_tmp.keys():
                        sparse_mtx = P_tmp.get(key)
                        P_tmp[key] = sparse.coo_matrix((sparse_mtx.data, (cells1_id[sparse_mtx.row], cells2_id[sparse_mtx.col])), shape=(n1,n2))
                        if not (key in P_cot.keys()):
                            P_cot[key]= P_tmp[key]
                        else:
                            P_cot[key].data = np.concatenate((P_cot[key].data,P_tmp[key].data))
                            P_cot[key].row  = np.concatenate((P_cot[key].row,P_tmp[key].row))
                            P_cot[key].col  = np.concatenate((P_cot[key].col,P_tmp[key].col))
                    
                    pbar.update()
    return P_cot

def sample_from_cluster(cell_cluster:pd.DataFrame,batchsize: int):
    '''
        Return the mini-batch samples according to cell_clusters.
        cell_cluster:pd.DataFrame
            cell_cluster['cluster']: numpy.ndarray
        Return: batches: list(list(int))
    '''
    clusters = [group.sample(frac=1).index.tolist() for _, group in cell_cluster.groupby('cluster')] # group and shuffle
    batch_nums = len(cell_cluster) // batchsize

    batches = [list() for _ in range(batch_nums)]
    for cluster in clusters:
        cells_per_batch = len(cluster) / batch_nums
        begin_loc = 0
        for i in range(batch_nums):
            end_loc = min(round((i+1)*cells_per_batch), len(cluster))
            batch_cells = cluster[begin_loc:end_loc]
            batches[i].extend(batch_cells)
            begin_loc = end_loc
    return batches


def transfer_anchor(S:pd.DataFrame, D:pd.DataFrame, A : pd.DataFrame, M: np.ndarray, \
                    cell1_cluster, cell2_cluster, pair = False, wnn = None, mnn = None, minibatch = 0, batchsize = 400,\
                    eps_p=1e-1, rho_mu=1e1, rho_nu=1e1, nitermax=1e3, stopthr=1e-8, device='cpu', q_cutoff = 0.75
                    )->pd.DataFrame:
    n1 = len(S.index)
    n2 = len(D.index)
    if (pair == True):
        if not (n1 == n2 and cell1_cluster.equals(cell2_cluster)):
            print("Please ensure that the data is pair!")
            return
        if wnn is None:
            print("Please input wnn similarity for paired data!")
            return
        # select anchors for paired data
        anchor_indnumsS = select_central_cells(wnn, list(cell1_cluster['cluster']), q_cutoff = q_cutoff)
        anchor_indnumsD = anchor_indnumsS
    else:
        if mnn is None:
            print("Please input mnn similarity for unpaired data!")
            return
        # select anchors for unpaired data
        anchor_indnumsS, anchor_indnumsD = select_anchor_cells_unpaired(S,D,mnn,list(cell1_cluster['cluster']),list(cell2_cluster['cluster']), q_cutoff = q_cutoff, scale = 0)    
    print(f'Identify {len(anchor_indnumsS),len(anchor_indnumsD)} anchor cells!')
    S_anchor_batch = S.iloc[anchor_indnumsS]
    D_anchor_batch = D.iloc[anchor_indnumsD]
    M_anchor_batch = M.iloc[anchor_indnumsS,anchor_indnumsD]
    cell1_cluster = cell1_cluster.iloc[anchor_indnumsS]
    cell2_cluster = cell2_cluster.iloc[anchor_indnumsD]
    cell1_cluster = cell1_cluster.reset_index(drop=True)
    cell2_cluster = cell2_cluster.reset_index(drop=True)
    r_sampleRate = len(anchor_indnumsS) / n1
    c_sampleRate = len(anchor_indnumsD) / n2
    rho_mu = rho_mu + eps_p * np.log(1/c_sampleRate)
    rho_nu = rho_nu + eps_p * np.log(1/r_sampleRate)
    P_anchor_tensor = fgot_sparse_tensor(S_anchor_batch, D_anchor_batch, A, M_anchor_batch,\
        cell1_cluster=cell1_cluster, cell2_cluster=cell2_cluster, minibatch=minibatch,\
        batchsize = batchsize, pair = pair, eps_p=eps_p, rho_mu=rho_mu, rho_nu=rho_nu,\
        nitermax=nitermax, stopthr=stopthr, device=device, fastMinibatch=True)
    P = fgot_tol(P_anchor_tensor)
    P = pd.DataFrame(P, index=S_anchor_batch.index, columns=D_anchor_batch.index)
    return P

def align(X, Y, P, mode = "ATAC2RNA"):
    # X: ATAC
    # Y: RNA
    # Projecting the X domain onto the Y domain
    if mode == "ATAC2RNA":
        Y1 = np.array(Y)
        weights = np.sum(P.T, axis = 0)
        X1 = np.matmul(P, np.array(Y))/ weights[:, None]
     # Projecting the Y domain onto the X domain
    elif mode == "RNA2ATAC":
        X1 = np.array(X)
        weights = np.sum(P, axis = 0)
        Y1 =np.matmul(P.T, np.array(X))/ weights[:, None]
    return X1, Y1


# for each gene, we extract cell type level regulatory intensity
# for each gene, average link intensity of the corresponding cells of the same group.
# golden standard: true link in feature matrix is 1, otherwise 0.
def fgot_analysis_link_intensity_for_each_celltype(P_tensor, feature_matrix, cellatac_cluster, cellrna_cluster, mode='ATAC2RNA'):
    # Transpose feature matrix if mode is RNA2ATAC
    if mode == 'RNA2ATAC':
        feature_matrix = feature_matrix.T
    # Build a dictionary to store gene pairs
    pairs = {}
    for key in tqdm(P_tensor.keys()):
        f_i, f_j = (key[0], key[1]) if mode == 'RNA2ATAC' else (key[1], key[0])
        pairs.setdefault(f_i, []).append(f_j)
    clusters = list(set(cellrna_cluster) & set(cellatac_cluster))
    intensity_s = []
    for f_i, f_j_list in tqdm(pairs.items()):
        int_i = np.zeros((len(clusters), len(f_j_list)))
        for i, cluster in enumerate(clusters):
            index1 = [k for k, x in enumerate(cellatac_cluster) if x == cluster]  # Index for ATAC
            index2 = [k for k, x in enumerate(cellrna_cluster) if x == cluster]  # Index for RNA
            for j, f_j in enumerate(f_j_list):
                key_j = (f_i, f_j) if mode == 'RNA2ATAC' else (f_j, f_i)
                sparse_matrix = P_tensor[key_j].todense()
                int_i[i, j] = sparse_matrix[np.ix_(index1, index2)].mean()
        # Convert to DataFrame and set columns as clusters (cell types)
        int_i_df = pd.DataFrame(int_i, index=clusters, columns=f_j_list)
        # Set row names as gene-peak combinations
        int_i_df = int_i_df.T
        int_i_df.index = [f"{f_i}-{peak}" for peak in int_i_df.index]
        intensity_s.append(int_i_df)
    # Combine all intensity DataFrames into one
    intensity_df = pd.concat(intensity_s, axis=0)
    return intensity_df


def refine(sample_id, pred, distance_mat, shape="hexagon"):
    '''
        Refine in single modality
        Thanks to https://github.com/jianhuupenn/SpaGCN/blob/master/SpaGCN_package/SpaGCN/util.py
        This function is inspired by SpaGCN.util.refine(sample_id,pred,dis,shape)
        For paired data with spatial information only.
    '''
    refined_pred=[]
    pred=pd.DataFrame({"pred": pred}, index=sample_id)
    distance_df=pd.DataFrame(distance_mat, index=sample_id, columns=sample_id) # 一个矩阵，每个元素是两个点之间的距离
    if shape=="hexagon":
        num_nbs=6
    elif shape=="square":
        num_nbs=4
    else:
        print("Shape not recongized, shape='hexagon' for Visium data, 'square' for ST data.")
    for i in range(len(sample_id)):
        index=sample_id[i]
        dis_tmp=distance_df.loc[index, :].sort_values()
        nbs=dis_tmp[0:num_nbs+1]  #choose the nearest neighbor points
        nbs_pred=pred.loc[nbs.index, "pred"]
        self_pred=pred.loc[index, "pred"]
        value_count=nbs_pred.value_counts()
        if (value_count.loc[self_pred]<num_nbs/2) and (np.max(value_count)>num_nbs/2):
            refined_pred.append(value_count.idxmax())
        else:           
            refined_pred.append(self_pred)
    return refined_pred
