import torch
import numpy as np
import pandas as pd
from scipy import sparse
from .utils import coo_submatrix_pull,fgot_flatten_matrix_sparse

def OTSolver(a:np.ndarray,
        b:np.ndarray,
        C:np.ndarray,
        eps_p:float,
        rho_mu:float,
        rho_nu:float,
        nitermax=1000,
        stopthr=1e-8,
        verbose=False,
        device='cpu'):
    """ The main function calling sinkhorn algorithm.

    Parameters
    ----------
    a : (ns,) numpy.ndarray
        Source distribution. The summation should be less than or equal to 1.
    b : (nt,) numpy.ndarray
        Target distribution. The summation should be less than or equal to 1.
    C : (ns,nt) numpy.ndarray
        The cost matrix possibly with infinity entries.
    eps_p :  float
        The coefficient of entropy regularization for P.
        = eps_mu : float, defaults to eps_p
        The coefficient of entropy regularization for mu.
        = eps_nu : float, defaults to eps_p
        The coefficient of entropy regularization for nu.
    rho mu : float
        The coefficient of penalty for unmatched mass for mu.
    rho nu : float
        The coefficient of penalty for unmatched mass for nu.
    nitermax : int, defaults to 10000
        The maximum number of iterations.
    stopthr : float, defaults to 1e-8
        The relative error threshold for stopping.
    verbose : boolean, defaults to False
        Whether to print algorithm logs.
    
    Edit records:
    Remove the option of 'solver', because we can only use sinkhorn(not 'momentum').
    Also, removed momentum parameters momentum_dt & momentum_beta.
    To use Sinkhorn algorithm, set eps_p=eps_mu=eps_nu.
    """
    # Return a zero matrix if either a or b is all zero
    nzind_a = np.where(a > 0)[0]
    nzind_b = np.where(b > 0)[0]
    if len(nzind_a) == 0 or len(nzind_b) == 0:
        P = sparse.coo_matrix(([],([],[])), shape=(len(a), len(b)))
        return P
    
    if device == 'cpu':
        P = OTSolver_sinkhorn_l1_sparse(a,b,C,eps_p,rho_mu,rho_nu, \
        nitermax=nitermax,stopthr=stopthr,verbose=verbose)
    else:
        P = OTSolver_sinkhorn_l1_sparse_torch(a,b,C,eps_p,rho_mu, rho_nu, \
            nitermax=nitermax,stopthr=stopthr,verbose=verbose,device=device)
    return P

def OTSolver_sinkhorn_l1_sparse(a,b,C,eps,rho_mu,rho_nu,nitermax=1000,stopthr=1e-6,verbose=False):
    """ Solve the unnormalized optimal transport with l1 penalty in sparse matrix format.

    Parameters
    ----------
    a : (ns,) numpy.ndarray
        Source distribution. The summation should be less than or equal to 1.
    b : (nt,) numpy.ndarray
        Target distribution. The summation should be less than or equal to 1.
    C : (ns,nt) scipy.sparse.coo_matrix
        The cost matrix in coo sparse format. The entries exceeds the cost cutoff are omitted. The naturally zero entries should be explicitely included.
    eps : float
        The coefficient for entropy regularization.
    rho_mu : float
        The coefficient for penalizing unmatched mass for mu.
    rho_nu : float
        The coefficient for penalizing unmatched mass for nu.
    nitermax : int, optional
        The max number of iterations. Defaults to 10000.
    stopthr : float, optional
        The threshold for terminating the iteration. Defaults to 1e-8.

    Returns
    -------
    (ns,nt) scipy.sparse.coo_matrix
        The optimal transport matrix. The locations of entries should agree with C and there might by explicit zero entries.
    """
    tmp_K = C.copy()
    f = np.zeros_like(a)
    g = np.zeros_like(b)
    r = np.zeros_like(a)
    s = np.zeros_like(b)
    niter = 0
    err = 100
    while niter <= nitermax and err > stopthr:
        fprev = f
        gprev = g
        # Iteration
        tmp_K.data = np.exp( ( -C.data + f[C.row] + g[C.col] ) / eps )
        f = eps * np.log(a) \
            - eps * np.log( np.sum( tmp_K, axis=1 ).A.reshape(-1) \
            + np.exp( ( -rho_mu + f ) / eps ) ) + f
        tmp_K.data = np.exp( ( -C.data + f[C.row] + g[C.col] ) / eps )
        g = eps * np.log(b) \
            - eps * np.log( np.sum( tmp_K, axis=0 ).A.reshape(-1) \
            + np.exp( ( -rho_nu + g ) / eps ) ) + g
        # Check relative error
        if niter % 10 == 0:
            err_f = abs(f - fprev).max() / max(abs(f).max(), abs(fprev).max(), 1.)
            err_g = abs(g - gprev).max() / max(abs(g).max(), abs(gprev).max(), 1.)
            err = 0.5 * (err_f + err_g)
        niter = niter + 1

    if verbose:
        print('Number of iterations in fgot:', niter)
    tmp_K.data = np.exp( ( -C.data + f[C.row] + g[C.col] ) / eps )
    return tmp_K


def OTSolver_sinkhorn_l1_sparse_torch(a, b, C, eps, rho_mu,rho_nu, nitermax=1000, stopthr=1e-6, verbose=False, device='cuda'):
    """ Solve the unnormalized optimal transport with l1 penalty in sparse matrix format using GPU.

    Parameters
    ----------
    a : (ns,) numpy.ndarray/torch.tensor
        Source distribution. The summation should be less than or equal to 1.
    b : (nt,) numpy.ndarray/torch.tensor
        Target distribution. The summation should be less than or equal to 1.
    C : (ns,nt) scipy.sparse.coo_matrix/torch.sparse_coo_tensor
        The cost matrix in coo sparse format. The entries exceeding the cost cutoff are omitted. The naturally zero entries should be explicitly included.
    eps : float
        The coefficient for entropy regularization.
    rho_mu : float
        The coefficient for penalizing unmatched mass for mu.
    rho_nu : float
        The coefficient for penalizing unmatched mass for nu.
    nitermax : int, optional
        The max number of iterations. Defaults to 10000.
    stopthr : float, optional
        The threshold for terminating the iteration. Defaults to 1e-8.

    Returns
    -------
    (ns,nt) scipy.sparse.coo_matrix
        The optimal transport matrix. The locations of entries should agree with C and there might be explicit zero entries.
    """
    # Convert the NumPy arrays and sparse matrix to PyTorch tensors and sparse tensor
    #print("Edited: Retrieve Cuda Memory in OTSolver_sinkhorn_torch")
    if device != 'cpu':
        torch.cuda.empty_cache()
    if(isinstance(a, np.ndarray)):
        a = torch.tensor(a, dtype=torch.float32, device=device)
    else:
        a = a.to(device)
    if(isinstance(b, np.ndarray)):
        b = torch.tensor(b, dtype=torch.float32, device=device)
    else:
        b = b.to(device)
    if(isinstance(C, sparse.coo_matrix)):
        C = torch.sparse_coo_tensor(torch.LongTensor(np.array([C.row, C.col])),
                                    torch.tensor(C.data, dtype=torch.float32),
                                    torch.Size(C.shape)).to(device)
    else:
        C = C.to(device)
    
    '''
    a : (ns,) torch.Tensor
    b : (nt,) torch.Tensor
    C : torch.sparse.FloatTensor
    '''
    
    tmp_K = C.clone()
    f = torch.zeros_like(a)
    g = torch.zeros_like(b)
    niter = 0
    err = 100

    while niter <= nitermax and err > stopthr:
        fprev = f.clone()
        gprev = g.clone()
        # Iteration
        tmp_K_value = torch.exp((-C._values() + f[C._indices()[0]] + g[C._indices()[1]]) / eps)
        tmp_K = torch.sparse_coo_tensor(C._indices(), tmp_K_value, C.shape)
        f = eps * torch.log(a) \
            - eps * torch.log(tmp_K.sum(dim=1).to_dense() + torch.exp((-rho_mu + f) / eps)) + f
        tmp_K_value = torch.exp((-C._values() + f[C._indices()[0]] + g[C._indices()[1]]) / eps)
        tmp_K = torch.sparse_coo_tensor(tmp_K._indices(), tmp_K_value, tmp_K.shape)
        g = eps * torch.log(b) \
            - eps * torch.log(tmp_K.sum(dim=0).to_dense() + torch.exp((-rho_nu + g) / eps)) + g
        #f[f == np.inf] = 0 # This would prolong the time very much
        # Sparse conditions do not produce an early end
        
        if niter % 10 == 0:
            err_f = torch.abs(f - fprev).max() / max(torch.abs(f).max(), torch.abs(fprev).max(), 1.)
            err_g = torch.abs(g - gprev).max() / max(torch.abs(g).max(), torch.abs(gprev).max(), 1.)
            err = 0.5 * (err_f + err_g)
        niter = niter + 1
    
    del err_f,err_g,err
    del fprev,gprev
    
    if verbose:
        print('Number of iterations in fgot:', niter)
    tmp_K_value = torch.exp((-C._values() + f[C._indices()[0]] + g[C._indices()[1]]) / eps)
    tmp_K = torch.sparse_coo_tensor(C._indices(), tmp_K_value, C.shape)
    del tmp_K_value,f,g
    
    indices = tmp_K._indices().cpu().numpy()
    values = tmp_K._values().cpu().numpy()
    
    if np.isinf(values).any():
        print("the solve of fgot_sinkhorn_l1_sparse_torch has inf, maybe please adjust the parameter eps_p to a higher value!")
    tmp_K = sparse.coo_matrix((values, (indices[0], indices[1])), shape=tmp_K.shape)
    del indices,values
    
    return tmp_K

def OTSolver_sinkhorn_l1_sparse_torch_step(a, b, C, f,g,eps, rho_mu,rho_nu, nitermax=500, stopthr=1e-6, verbose=False, device='cuda'):
    """ Step to unnormalized optimal transport with l1 penalty in sparse matrix format using GPU.

    Parameters
    ----------
    a : (ns,) numpy.ndarray/torch.tensor
        Source distribution. The summation should be less than or equal to 1.
    b : (nt,) numpy.ndarray/torch.tensor
        Target distribution. The summation should be less than or equal to 1.
    C : (ns,nt) scipy.sparse.coo_matrix/torch.sparse_coo_tensor
        The cost matrix in coo sparse format. The entries exceeding the cost cutoff are omitted. The naturally zero entries should be explicitly included.
    eps : float
        The coefficient for entropy regularization.
    rho_mu : float
        The coefficient for penalizing unmatched mass for mu.
    rho_nu : float
        The coefficient for penalizing unmatched mass for nu.
    nitermax : int, optional
        The max number of iterations. Defaults to 10000.
    stopthr : float, optional
        The threshold for terminating the iteration. Defaults to 1e-8.

    Returns
    -------
    (ns,nt) scipy.sparse.coo_matrix
        The optimal transport matrix. The locations of entries should agree with C and there might be explicit zero entries.
    """
    # Convert the NumPy arrays and sparse matrix to PyTorch tensors and sparse tensor
    #print("Edited: Retrieve Cuda Memory in OTSolver_sinkhorn_torch")
    if device != 'cpu':
        torch.cuda.empty_cache()
    if(isinstance(a, np.ndarray)):
        a = torch.tensor(a, dtype=torch.float32, device=device)
    else:
        a = a.to(device)
    if(isinstance(b, np.ndarray)):
        b = torch.tensor(b, dtype=torch.float32, device=device)
    else:
        b = b.to(device)
    if(isinstance(C, sparse.coo_matrix)):
        C = torch.sparse_coo_tensor(torch.LongTensor(np.array([C.row, C.col])),
                                    torch.tensor(C.data, dtype=torch.float32),
                                    torch.Size(C.shape)).to(device)
    else:
        C = C.to(device)
    if(isinstance(f, np.ndarray)):
        f = torch.tensor(f, dtype=torch.float32, device=device)
    else:
        f = f.to(device)
    if(isinstance(g, np.ndarray)):
        g = torch.tensor(g, dtype=torch.float32, device=device)
    else:
        g = g.to(device)
    
    '''
    a : (ns,) torch.Tensor
    b : (nt,) torch.Tensor
    C : torch.sparse.FloatTensor
    '''
    niter = 0
    err = 100
    while niter <= nitermax and err > stopthr:
        fprev = f.clone()
        gprev = g.clone()
        tmp_K = C.clone()
        # Iteration
        tmp_K_value = torch.exp((-C._values() + f[C._indices()[0]] + g[C._indices()[1]]) / eps)
        tmp_K = torch.sparse_coo_tensor(C._indices(), tmp_K_value, C.shape)
        f = eps * torch.log(a) \
            - eps * torch.log(tmp_K.sum(dim=1).to_dense() + torch.exp((-rho_mu + f) / eps)) + f
        tmp_K_value = torch.exp((-C._values() + f[C._indices()[0]] + g[C._indices()[1]]) / eps)
        tmp_K = torch.sparse_coo_tensor(tmp_K._indices(), tmp_K_value, tmp_K.shape)
        g = eps * torch.log(b) \
            - eps * torch.log(tmp_K.sum(dim=0).to_dense() + torch.exp((-rho_nu + g) / eps)) + g
        
        #f[f == np.inf] = 0 # This would prolong the time very much
        # Sparse conditions do not produce an early end
        
        tmp_K_value = torch.exp((-C._values() + f[C._indices()[0]] + g[C._indices()[1]]) / eps)
        tmp_K = torch.sparse_coo_tensor(C._indices(), tmp_K_value, C.shape)
        del tmp_K_value
        if niter % 10 == 0:
            err_f = torch.abs(f - fprev).max() / max(torch.abs(f).max(), torch.abs(fprev).max(), 1.)
            err_g = torch.abs(g - gprev).max() / max(torch.abs(g).max(), torch.abs(gprev).max(), 1.)
            err = 0.5 * (err_f + err_g)
        niter = niter + 1
    del err_f,err_g
    del fprev,gprev
    
    f_ret = f.cpu().numpy()
    g_ret = g.cpu().numpy()
    err_ret = err.cpu().numpy()
    del f,g,err
    
    indices = tmp_K._indices().cpu().numpy()
    values = tmp_K._values().cpu().numpy()
    tmp_K = sparse.coo_matrix((values, (indices[0], indices[1])), shape=tmp_K.shape)
    del indices,values
    
    return tmp_K,f_ret,g_ret,err_ret


def fgot_flatten_solve(S:pd.DataFrame, D:pd.DataFrame, A:pd.DataFrame, M:pd.DataFrame,\
    eps_p=1e-1, rho_mu=1e1, rho_nu=1e1 , nitermax=1e3, stopthr=1e-8, device='cpu',\
    c_sampleRate=1,r_sampleRate=1)->dict:
    """ Solve the feature-guided optimal transport problem in sparse format.
    Translate the problem to 2D OT problem and solve it with FGOT solver.
    
    Parameters
    ----------
    S : (n_pos_s,ns_s) panda.DataFrame
        Source distributions over `n_pos_s` positions of `ns_s` source species.
    D : (n_pos_d,ns_d) panda.DataFrame
        Destination distributions over `n_pos_d` positions of `ns_d` destination species.
    A : (ns_s,ns_d) panda.DataFrame
        The cost coefficients for source-destination species pairs. An infinity value indicates that the two species cannot be coupled.
    M : (n_pos_s,n_pos_d) panda.DataFrame
        The distance (cost) matrix among the positions.
    eps_p : float, defaults to 1e-1
        The coefficient for entropy regularization of P.
    rho : float, defaults to 1e2
        The coefficient for penalizing unmatched mass.
    nitermax : int, optional
        The maximum number of iterations in the unormalized OT problem. Defaults to 1e3.
    stopthr : float, optional
        The relatitive error threshold for terminating the iteration. Defaults to 1e-8.
    device  : str, optional
        To run the algorithm on which device. Defaults to 'cpu'.
        Choose from 'cpu' and 'cuda'/'cuda:x'.
    adjust  : bool, optional
        Whether to adjust the rho_mu, rho_nu, a and b according to the sample rate. Defaults to True.
    c_sampleRate : float, optional
        The sample rate of the destination species. Defaults to 1.
    r_sampleRate : float, optional
        The sample rate of the source species. Defaults to 1.
    
    Returns
    -------
    A dictionary of scipy.sparse.coo_matrix
        The transport plan in coo sparse format for source species i and destinaton species j can be retrieved with the key (i,j).
    """
    gene_names = list(S.columns)
    peak_names = list(D.columns)
    S = np.array(S)
    D = np.array(D)
    A = np.array(A)
    M = np.array(M)
    
    np.set_printoptions(precision=2)
    n_pos_s, ns_s = S.shape
    n_pos_d, ns_d = D.shape
    S[S<0] = 0
    D[D<0] = 0
    max_amount = max( S.sum(), D.sum() )
    S = S / max_amount
    D = D / max_amount
    
    # Set up the large collective OT problem
    ##############################################################################
    # flatten sparse
    C = fgot_flatten_matrix_sparse(M, S, D, A, ns_s, ns_d, n_pos_s, n_pos_d)
    a = S.flatten('F')
    b = D.flatten('F')
    ##################################################################################################
    
    
    # Solve the problem on nonzero mass
    nzind_a = np.where(a > 0)[0]
    nzind_b = np.where(b > 0)[0]
    C_nz = coo_submatrix_pull(C, nzind_a, nzind_b)
    del C

    a = a*c_sampleRate
    b = b*r_sampleRate
    tmp_P = OTSolver(a[nzind_a], b[nzind_b], C_nz, eps_p, rho_mu, rho_nu,\
        nitermax=nitermax, stopthr=stopthr, device=device)
    del C_nz

    P = sparse.coo_matrix((tmp_P.data, (nzind_a[tmp_P.row], nzind_b[tmp_P.col])), shape=(len(a),len(b)))
    P = P.tocsr()

    # Output a dictionary of transport plans & check the result
    P_expand = {}
    for i in range(ns_s):
        for j in range(ns_d):
            if not np.isinf(A[i,j]):
                tmp_P = P[i*n_pos_s:(i+1)*n_pos_s, j*n_pos_d:(j+1)*n_pos_d]
                P_expand[(gene_names[i],peak_names[j])] = tmp_P.tocoo() * max_amount
                sparse_mtx = P_expand.get((gene_names[i],peak_names[j]))
                if np.isnan(sparse_mtx.data).any():
                    print("Error! The model parameters are inappropriate and please consider enlarging the value of the parameters esp_p!")
                    return
                if np.any(sparse_mtx.data<0):
                    print("the solution has number < 0!")
                    return
    
    return P_expand