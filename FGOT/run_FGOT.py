import argparse
import numpy as np
import pandas as pd
import matplotlib.pylab as pl
import scanpy as sc
from anndata import AnnData
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import FGOT
from FGOT import preprocess as pre
from FGOT.fgot import fgot_sparse_tensor, fgot_tol, align
import pickle

def parse_args():
    parser = argparse.ArgumentParser(description="Run FGOT analysis")
    parser.add_argument("-i", "--input_type", type=str, required=True, help="Input type (e.g., paired or unpaired)")
    parser.add_argument("--rna_data", type=str, required=True, help="Path to the normalized RNA data file")
    parser.add_argument("--atac_data", type=str, required=True, help="Path to the normalized ATAC data file")
    parser.add_argument("--cluster_info", type=str, required=True, help="Path to the RNA cluster info file")
    parser.add_argument("--atac_cluster_info", type=str, default=None, help="Path to the ATAC cluster info file (optional, required if unpaired)")
    parser.add_argument("--feature_matrix", type=str, required=True, help="Path to the peak-gene matrix file (where a value of 1 indicates that peaks are located close to the gene)")
    parser.add_argument("--corr_matrix", type=str, required=True, help="Path to the cell corr matrix file (snn or wnn matrix)")
    parser.add_argument("--output_prefix", type=str, required=True, help="Prefix for output files")
    parser.add_argument("--device", type=str, default="cuda:0", help="Device to use (e.g., cuda:0 or cpu)")
    parser.add_argument("--minibatch", type=int, default=0, help="Whether to use minibatch (default 0 for not use)")
    parser.add_argument("--batchsize", type=int, default=400, help="Batch size for mini-batch computing (optional, required if minibatch = 1)")
    parser.add_argument("--eps_p", type=float, default=1e-2, help="Epsilon for tensor computation")
    parser.add_argument("--rho_mu", type=float, default=1e1, help="Parameter rho_mu")
    parser.add_argument("--rho_nu", type=float, default=1e1, help="Parameter rho_nu")
    parser.add_argument("--nitermax", type=int, default=1000, help="Maximum number of iterations")
    parser.add_argument("--stopthr", type=float, default=1e-8, help="Stopping threshold")
    parser.add_argument("--fastMinibatch", action='store_false', default=True, help="Enable fast mini-batch mode if default True, else anchor-based mini-batch strategy (used when minibatch = 1)")
    parser.add_argument("--lam", type=float, default=0.5, help="Lambda for regularization")
    parser.add_argument("--adjust", action='store_false', default=True, help="Whether to adjust the regularization parameters based on the sampling rates in mini-batch mode (default True)")
    parser.add_argument("--align_mode", type=str, default="ATAC2RNA", help="alignment mode (default \"ATAC2RNA\", else \"RNA2ATAC\")")
    return parser.parse_args()


def main():
    args = parse_args()

    # Load data
    RNA_data = pd.read_csv(args.rna_data, sep=',', index_col=0)
    ATAC_data = pd.read_csv(args.atac_data, sep=',', index_col=0)
    RNA_cluster = pd.read_csv(args.cluster_info, sep='\t')
    RNA_cluster = RNA_cluster.reset_index(drop=True)
    RNA_cluster = RNA_cluster.rename(columns={"cell_id": "cell", "celltype":"cluster"})
    if args.input_type == "paired":
        ATAC_cluster = RNA_cluster.copy()
    else:
        if args.atac_cluster_info is not None:
            ATAC_cluster = pd.read_csv(args.atac_cluster_info, sep="\t")
            ATAC_cluster = ATAC_cluster.reset_index(drop=True)
            ATAC_cluster = ATAC_cluster.rename(columns={"cell_id": "cell", "celltype":"cluster"})
        else:
            raise ValueError("ATAC cluster info file is required when input type is unpaired.")
    X1 = ATAC_data.T
    X2 = RNA_data.T
    (n1, d1) = X1.shape
    (n2, d2) = X2.shape
    print("Dimensions of input datasets are: ", "X1 = ", X1.shape, " X2 = ", X2.shape)
    basic_labels1 = np.array(ATAC_cluster['cluster'])
    basic_labels2 = np.array(RNA_cluster['cluster'])
    print("Dimensions of basic_labels are: ", "basic_labels1 = ", basic_labels1.shape, " basic_labels2 = ", basic_labels2.shape)
    cell_names1 = list(ATAC_data.columns)
    cell_names2 = list(RNA_data.columns)
    peak_names = list(ATAC_data.index)
    gene_names = list(RNA_data.index)

    # build a priori feature graph
    print(f"Building the priori feature graph for FGOT...")
    # for simulation multi-omics data
    feature_matrix = pd.read_csv(args.feature_matrix, sep='\t', index_col=0)
    feature_matrix = feature_matrix.replace(0, np.inf)
    # for true multi-omics data which we have gene annotations（4 columns: chr, starts, ends, genes）
    # promoters = pd.read_csv('hg38.promoter.regions.txt', sep = '\t') # hg19/mm10
    # feature_matrix = pre.prior_feature_graph(promoters, peak_names, gene_names)
    print("Feature matrix loaded with shape: ", feature_matrix.shape)
    
    # compute cost across modalities
    print(f"Computing the cost matrix across modalities for FGOT...")
    if args.input_type == "paired":
        wnn = pd.read_csv(args.corr_matrix, sep='\t', index_col=0)
        cost = np.array(np.exp(1 - wnn))
    else:
        S_m = pd.read_table(args.corr_matrix, index_col=0).T
        S_sm = pre.smooth_cell_similarity_byLaplacian1(S_m, X1, X2)
        S = S_sm/np.max(S_sm) + S_m/np.max(S_m)
        cost = np.array(np.exp(np.max(S) - S))
    cost = cost - np.min(cost)
    cost = pd.DataFrame(cost, index = cell_names1, columns=cell_names2)
    print("Cost matrix loaded with shape: ", cost.shape)
    
    # optional scale
    # If you want more accurate regulatory results for each gene, avoid scaling. 
    # If you prefer to save memory and get faster results for large-scale data, consider scaling.
    
    # print(f"Scaling the input data...")
    # scaler = StandardScaler()
    # X1, X2 = scaler.fit_transform(X1), scaler.fit_transform(X2)
    # X1 = pd.DataFrame(X1,index= cell_names1,columns=peak_names)
    # X2 = pd.DataFrame(X2,index= cell_names2,columns=gene_names)
    
    # solve the feature-guided optimal transport(FGOT) problem
    print(f"Solve the feature-guided optimal transport(FGOT) problem...")
    if args.input_type == "paired":
        pair = True
    else:
        pair = False
    P_tensor = fgot_sparse_tensor(X1, X2, feature_matrix, cost,\
        ATAC_cluster, RNA_cluster, minibatch=args.minibatch, batchsize=args.batchsize, \
        pair=pair, eps_p=args.eps_p, rho_mu=args.rho_mu, rho_nu=args.rho_nu, \
        nitermax=args.nitermax, stopthr=args.stopthr, device=args.device, \
        fastMinibatch=args.fastMinibatch, lam=args.lam, adjust=args.adjust)
    pickle.dump(P_tensor, open(f"{args.output_prefix}P_tensor.pickle", "wb"))
    print(f"FGOT tensor saved as {args.output_prefix}P_tensor.pickle")
  
    # make alignment
    print(f"Making alignment...")
    P = fgot_tol(P_tensor)
    X1_aligned, X2_aligned = align(X1, X2, P, args.align_mode)
    data_aligned = np.concatenate((X2_aligned, X1_aligned), axis=0)
    adata_aligned = AnnData(data_aligned)
    adata_aligned.obs['batch'] = np.array(['RNA'] * n2 + ['ATAC'] * n1)
    adata_aligned.obs['label'] = np.concatenate((basic_labels2, basic_labels1),axis=0)
    sc.tl.pca(adata_aligned)
    sc.pp.neighbors(adata_aligned, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata_aligned)
    sc.pl.umap(adata_aligned, color=['batch', 'label'],s = 80, show=False)
    pl.savefig(f"{args.output_prefix}alignment_UMAP.pdf")
    print(f"UMAP figure saved as {args.output_prefix}alignment_UMAP.pdf")
    
    ## infer the celltype specific regulatory intensity from the transport tensor
    print(f"Analysising celltype specific regulatory intensity...")
    intensity_df = FGOT.fgot.fgot_analysis_link_intensity_for_each_celltype(P_tensor, feature_matrix, basic_labels1, basic_labels2, mode = args.align_mode)
    intensity_df.to_csv(f"{args.output_prefix}celltype_specific_regulatory_df.csv", index=True, sep=",")
    print(f"Regulatory intensity dataframe saved as {args.output_prefix}celltype_specific_regulatory_df.csv")

    print("FGOT analysis complete.")

if __name__ == "__main__":
    main()


# shell command
# python3 /home/nas2/biod/yangchenghui/FGOT-master/FGOT/run_FGOT.py -i paired \
# --rna_data /home/nas2/biod/yangchenghui/FGOT_scMultiSim_simulation/simulation_data_discrete_4types/normalized_rna_data.txt \
# --atac_data /home/nas2/biod/yangchenghui/FGOT_scMultiSim_simulation/simulation_data_discrete_4types/normalized_atac_data.txt \
# --cluster_info /home/nas2/biod/yangchenghui/FGOT_scMultiSim_simulation/simulation_data_discrete_4types/celltype_info.txt \
# --feature_matrix /home/nas2/biod/yangchenghui/FGOT_scMultiSim_simulation/simulation_data_discrete_4types/feature_matrix.txt \
# --corr_matrix /home/nas2/biod/yangchenghui/FGOT_scMultiSim_simulation/simulation_data_discrete_4types/simu_wnn.txt \
# --output_prefix /home/nas2/biod/yangchenghui/FGOT_scMultiSim_simulation/FGOT_output_test/
