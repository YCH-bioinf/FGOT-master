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
from .utils import select_confident_cells,compute_D_kl


def fgot_tol(P_4dim):
    shape = P_4dim.get(list(P_4dim.keys())[0]).shape
    P_2dim = np.zeros(shape)
    for key in tqdm(P_4dim.keys()):
        sparse_matrix = P_4dim[key].todense()
        P_2dim += sparse_matrix
    return P_2dim


def fgot_sparse_4dim(S:pd.DataFrame, D:pd.DataFrame, A : pd.DataFrame, M: np.ndarray, \
                    cell1_cluster=None, cell2_cluster=None, minibatch=0, batchsize=100, pair = False,\
                    eps_p=1e-1, rho=1e1, nitermax=1e3, stopthr=1e-8, device='cpu',\
                    fastMinibatch=True, P_anchor:pd.DataFrame=None, lam:float=0.01,adjust=True)->dict:
    """ 
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
            eps_p=eps_p, rho_mu=rho, rho_nu=rho,\
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
            from .utils import get_snn_mat
            Snn1 = get_snn_mat(S)
            Snn2 = get_snn_mat(D)
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
                    rho_mu = rho + eps_p * np.log(1/c_sampleRate)
                    rho_nu = rho + eps_p * np.log(1/r_sampleRate)
                else:
                    r_sampleRate = 1
                    c_sampleRate = 1
                    rho_mu = rho
                    rho_nu = rho
                
                if fastMinibatch==False:
                    Snn1_batch = Snn1[cells1_id,:][:,anchors1_id]
                    Snn2_batch = Snn2[cells2_id,:][:,anchors2_id]
                    for i in range(len(cells1_id)):
                        for j in range(len(cells2_id)):
                            M_batch.iloc[i,j] = M_batch.iloc[i,j] + lam * np.sum(np.multiply(compute_D_kl(Snn1_batch, Snn2_batch, i, j),P_anchor))
                
                P_tmp = fgot_flatten_solve(S_batch, D_batch, A, M_batch,\
                    eps_p=eps_p, rho_mu=rho_mu, rho_nu=rho_nu,\
                    nitermax=nitermax, stopthr=stopthr,device=device,\
                    adjust=adjust,c_sampleRate=c_sampleRate,r_sampleRate=r_sampleRate)
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
                        rho_mu = rho + eps_p * np.log(1/c_sampleRate)
                        rho_nu = rho + eps_p * np.log(1/r_sampleRate)
                    else:
                        r_sampleRate = 1
                        c_sampleRate = 1
                        rho_mu = rho
                        rho_nu = rho
                    
                    if fastMinibatch==False:
                        Snn1_batch = Snn1[cells1_id,:][:,anchors1_id]
                        Snn2_batch = Snn2[cells2_id,:][:,anchors2_id]
                        for i in range(len(cells1_id)):
                            for j in range(len(cells2_id)):
                                M_batch.iloc[i,j] = M_batch.iloc[i,j] + lam * np.sum(np.multiply(compute_D_kl(Snn1_batch, Snn2_batch, i, j),P_anchor))
                    
                    P_tmp = fgot_flatten_solve(S_batch, D_batch, A, M_batch,\
                        eps_p=eps_p, rho_mu=rho_mu, rho_nu=rho_nu,\
                        nitermax=nitermax, stopthr=stopthr,device=device,\
                        adjust=adjust,c_sampleRate=c_sampleRate,r_sampleRate=r_sampleRate)
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

def adjust(X1_aligned:np.ndarray,X2_aligned:np.ndarray,paired=False,k=10):
    '''
    Adjust the X2_aligned to X1_aligned by translation and scaling.
    X1_aligned: np.ndarray
        The aligned source data.
    X2_aligned: np.ndarray
        The aligned target data.
    paired: bool, optional
        Whether the two datasets are paired. Defaults to False.
    k: int, optional
        The number of nearest neighbors for each randomly chosen cell. Defaults to 10.
    '''
    if paired:
        X1_mnn = X1_aligned
        X2_mnn = X2_aligned
    else:
        nn_X1 = NearestNeighbors(n_neighbors=k).fit(X1_aligned)
        nn_X2 = NearestNeighbors(n_neighbors=k).fit(X2_aligned)
        distances_X1_to_X2, indices_X1_to_X2 = nn_X1.kneighbors(X2_aligned)
        distances_X2_to_X1, indices_X2_to_X1 = nn_X2.kneighbors(X1_aligned)

        mnn_pairs_X1 = []
        mnn_pairs_X2 = []
        for id_x2 in range(indices_X1_to_X2.shape[0]):
            list_x1 = indices_X1_to_X2[id_x2]
            for id_x1 in list_x1:
                if id_x2 in indices_X2_to_X1[id_x1]:
                    if id_x1 not in mnn_pairs_X1:
                        mnn_pairs_X1.append(id_x1)
                    if id_x2 not in mnn_pairs_X2:
                        mnn_pairs_X2.append(id_x2)
        X1_mnn = X1_aligned[mnn_pairs_X1]
        X2_mnn = X2_aligned[mnn_pairs_X2]

    scale_mnn = np.abs((np.max(X1_mnn, axis=0) - np.min(X1_mnn, axis=0)) / (np.max(X2_mnn, axis=0) - np.min(X2_mnn, axis=0)))
    scale_mnn = np.where(np.isfinite(scale_mnn), scale_mnn, 1)
    X2_scaled_mnn = X2_mnn * scale_mnn
    mean_X2_mnn = np.mean(X2_scaled_mnn, axis=0)
    mean_X1_mnn = np.mean(X1_mnn, axis=0)
    translation_mnn = mean_X1_mnn - mean_X2_mnn
    X2_adjust = X2_aligned* scale_mnn + translation_mnn
    return X2_adjust, scale_mnn, translation_mnn

def cross_refine(sample_id, pred, distance_mat, shape="hexagon"):
    '''
        /*Warning： This function is still under test.*/
        Refine for two modalities with spatial information.
        Thanks to https://github.com/jianhuupenn/SpaGCN/blob/master/SpaGCN_package/SpaGCN/util.py
        This function is inspired by SpaGCN.util.refine(sample_id,pred,dis,shape)
        For paired data with spatial information only.
    '''
    refined_pred=[]
    pred1 = pred[:len(sample_id)]
    pred2 = pred[len(sample_id):]
    pred=pd.DataFrame({"pred_modal1": pred1,"pred_modal2":pred2}, index=sample_id)
    distance_df=pd.DataFrame(distance_mat, index=sample_id, columns=sample_id)
    if shape=="hexagon":
        num_nbs=6
    elif shape=="square":
        num_nbs=4
    else:
        print("Shape not recongized, shape='hexagon' for Visium data, 'square' for ST data.")
    for i in range(len(sample_id)):
        index=sample_id[i]
        dis_tmp=distance_df.loc[index, :].sort_values()
        nbs=dis_tmp[0:num_nbs+1] #choose the nearest neighbor points
        nbs_pred1=pred.loc[nbs.index, "pred_modal1"]
        nbs_pred2=pred.loc[nbs.index, "pred_modal2"]
        self_pred=pred.loc[index, "pred_modal1"]
        value_count=nbs_pred1.value_counts()+nbs_pred2.value_counts()
        if (value_count.loc[self_pred]<num_nbs/2) and (np.max(value_count)>num_nbs*2/2):
            refined_pred.append(value_count.idxmax())
        else:           
            refined_pred.append(self_pred)
    refined_pred = refined_pred + refined_pred
    return refined_pred

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

def transfer_anchor(S:pd.DataFrame, D:pd.DataFrame, A : pd.DataFrame, M: np.ndarray, \
                     cell1_cluster, cell2_cluster,\
                     eps_p=1e-1, rho=1e1, nitermax=1e3, stopthr=1e-8, device='cpu'
                     )->pd.DataFrame:
    n1 = len(S.index)
    n2 = len(D.index)
    ## Compute P anchor
    adataS= AnnData(S)
    adataS.obs['cell_type'] = list(cell1_cluster['cluster'])
    anchor_namesS,anchor_indnumsS = select_confident_cells(adataS,'cell_type')
    del adataS
    adataD= AnnData(D)
    adataD.obs['cell_type'] = list(cell2_cluster['cluster'])
    anchor_namesD,anchor_indnumsD = select_confident_cells(adataD,'cell_type')
    del adataD
    
    print(f'Identify {len(anchor_indnumsS),len(anchor_indnumsD)} anchor cells!')
    S_anchor_batch = S.iloc[anchor_indnumsS]
    D_anchor_batch = D.iloc[anchor_indnumsD]
    M_anchor_batch = M.iloc[anchor_indnumsS,anchor_indnumsD]
    r_sampleRate = len(anchor_indnumsS) / n1
    c_sampleRate = len(anchor_indnumsD) / n2
    rho_mu = rho + eps_p * np.log(1/c_sampleRate)
    rho_nu = rho + eps_p * np.log(1/r_sampleRate)
    P_anchor_4dim = fgot_flatten_solve(S_anchor_batch, D_anchor_batch, A, M_anchor_batch,\
        eps_p=eps_p, rho_mu=rho_mu, rho_nu=rho_nu,\
        nitermax=nitermax, stopthr=stopthr,device=device)
    P = fgot_tol(P_anchor_4dim)
    P = pd.DataFrame(P,index=anchor_namesS, columns=anchor_namesD)
    return P

def align(Y, X, P, X2Y=True, ajust = False):
    # Projecting the X domain onto the Y domain
    if X2Y:
        Y1 = np.array(Y)
        weights = np.sum(P, axis = 0)
        X1 = np.matmul(P.T, np.array(Y))/ weights[:, None]
        if ajust:
            X1,_,_ = adjust(Y1, X1)
    else:
        X1 = np.array(X)
        weights = np.sum(P.T, axis = 0)
        Y1 =np.matmul(P, np.array(X))/ weights[:, None]
        if ajust:
            Y1,_,_ = adjust(X1,Y1)
    return Y1, X1