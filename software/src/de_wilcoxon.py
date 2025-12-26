import argparse
import pandas as pd
import polars as pl
import numpy as np
import scanpy as sc
import time
from anndata import AnnData
from scipy.sparse import csr_matrix

np.random.seed(0)

def log_message(message, status="INFO"):
    """Logs messages in a structured format."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [{status}] {message}")

def main(expr_path, meta_path, group1, group2, output, condition_col, 
    padj_cutoff, logfc_cutoff, error_output):
    
    log_message("Starting Differential Expression analysis", "STEP")

    # 1. LOAD METADATA (Pandas is preferred for metadata as per repo style)
    df_meta = pd.read_csv(meta_path)

    # First check if we have enough data to run DE. Return errors otherwise
    errorLogs = []
    for group_ in [group1, group2]:
        sumGroup = sum(df_meta[condition_col] == group_)
        if sumGroup == 0:
            errorL = f"Statistics could not be calculated for group(s) {group_} because they are not present in any sample."
            print(errorL)
            errorLogs += [errorL]
    if len(errorLogs) > 0:
        # output empty tables
        cols = ["Contrast", "gene", "logfc", "pval", "pval_adj", "score", "minlog10padj", "regulationDirection"]
        df_empty = pd.DataFrame(columns=cols)
        df_empty.to_csv(output, index=False)
        df_empty.to_csv(output.replace('.csv', '_filtered.csv'), index=False)

        # output error logs table
        df_error = pd.DataFrame(columns=["Error", "value"])
        df_error['Error'] = errorLogs
        df_error['value'] = errorLogs
        df_error.to_csv(error_output, index=False)
        return
    else:
        df_error = pd.DataFrame(columns=["Error", "value"])
        df_error.to_csv(error_output, index=False)

    # 2. LOAD EXPRESSION DATA WITH POLARS & CATEGORICAL OPTIMIZATION
    log_message("Loading expression data with Polars optimization", "STEP")
    
    # Peek schema for flexible headers and identify columns for Categorical casting
    temp_scan = pl.scan_csv(expr_path)
    file_schema = temp_scan.collect_schema()
    column_names = set(file_schema.keys())

    schema_overrides = {
        "Sample": pl.Categorical,
        "Ensembl Id": pl.Categorical,
        "Cell Barcode": pl.Categorical,
        "Cell ID": pl.Categorical,
    }
    # Only apply overrides for columns that actually exist in the file
    schema_overrides = {k: v for k, v in schema_overrides.items() if k in column_names}
    
    df_expr_pl = pl.read_csv(expr_path, schema_overrides=schema_overrides)

    # Normalize header labels to support migration from 'Cell Barcode' to 'Cell ID'
    if 'Cell Barcode' not in df_expr_pl.columns and 'Cell ID' in df_expr_pl.columns:
        df_expr_pl = df_expr_pl.rename({'Cell ID': 'Cell Barcode'})

    # Create a unique identifier for each cell (using established SEPARATOR)
    SEPARATOR = '|||'
    df_expr_pl = df_expr_pl.with_columns(
        (pl.col('Sample').cast(str) + pl.lit(SEPARATOR) + pl.col('Cell Barcode').cast(str))
        .cast(pl.Categorical)
        .alias('UniqueCellId')
    )

    # 3. CONSTRUCT SPARSE MATRIX DIRECTLY (REFINED STYLE)
    log_message("Creating sparse matrix from long format data", "STEP")
    
    # Extract integer codes directly from categorical columns (much faster than join)
    row_codes_raw = df_expr_pl['UniqueCellId'].to_physical().to_numpy()
    col_codes_raw = df_expr_pl['Ensembl Id'].to_physical().to_numpy()
    # Use float32 for expression values (Scanpy standard)
    expression_values = df_expr_pl['Raw gene expression'].cast(pl.Float32).to_numpy()

    # Remap codes to 0-indexed contiguous using np.unique (efficient integer-based mapping)
    u_row_phys, row_idx = np.unique(row_codes_raw, return_inverse=True)
    u_col_phys, col_idx = np.unique(col_codes_raw, return_inverse=True)

    # Map labels to sorted ranks and get sorted unique IDs (efficiently processing only unique labels)
    unique_cell_ids, row_map = np.unique(df_expr_pl['UniqueCellId'].cat.get_categories().gather(u_row_phys).to_numpy(), return_inverse=True)
    unique_gene_ids, col_map = np.unique(df_expr_pl['Ensembl Id'].cat.get_categories().gather(u_col_phys).to_numpy(), return_inverse=True)

    # Final row and column codes are the mapped indices
    row_codes = row_map[row_idx].astype(np.int32)
    col_codes = col_map[col_idx].astype(np.int32)

    # Delete Polars objects to free memory
    del df_expr_pl, row_codes_raw, col_codes_raw, u_row_phys, u_col_phys, row_idx, col_idx, row_map, col_map
    
    # Pre-populate obs with Sample and Cell Barcode vectorially for efficient processing
    obs_df = pd.DataFrame(index=unique_cell_ids)
    split_ids = pd.Series(unique_cell_ids).str.split(SEPARATOR, n=1, expand=True, regex=False)
    obs_df['Sample'] = split_ids[0].values
    obs_df['Cell Barcode'] = split_ids[1].values

    # Create the sparse matrix and AnnData object
    adata = AnnData(
        X=csr_matrix((expression_values, (row_codes, col_codes)), shape=(len(unique_cell_ids), len(unique_gene_ids)), dtype=np.float32),
        obs=obs_df,
        var=pd.DataFrame(index=unique_gene_ids)
    )
    log_message(f"AnnData object created: {adata.n_obs} cells x {adata.n_vars} genes", "DONE")

    # 4. ALIGN METADATA
    log_message("Aligning metadata to cells", "STEP")
    ## Cases where we do DE using sample metadata
    if 'Cell ID' not in df_meta.columns:
        # Join on Sample while preserving the original index and row order
        adata.obs = adata.obs.join(df_meta.set_index('Sample'), on='Sample')
    ## Cases where we do DE using cell metadata
    else:
        meta_temp = df_meta.rename(columns={'Cell ID': 'Cell Barcode'})
        # Create matching UniqueCellId for merging
        meta_temp['UniqueCellId'] = meta_temp['Sample'].astype(str) + SEPARATOR + meta_temp['Cell Barcode'].astype(str)
        # Join on UniqueCellId while preserving the original index and row order
        # Drop redundant sample/barcode from metadata to avoid conflicts
        adata.obs = adata.obs.join(
            meta_temp.set_index('UniqueCellId').drop(columns=['Sample', 'Cell Barcode'], errors='ignore')
        )
    
    # 5. PREPROCESSING & ANALYSIS
    log_message("Normalizing and log-transforming", "STEP")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Subset to relevant groups
    adata.obs[condition_col] = adata.obs[condition_col].astype(str)
    adata = adata[adata.obs[condition_col].isin([group1, group2])].copy()

    # Run DE test
    log_message(f"Running Wilcoxon rank-sum test: {group1} vs {group2}", "STEP")
    sc.tl.rank_genes_groups(
        adata,
        groupby=condition_col,
        groups=[group1],
        reference=group2,
        method='wilcoxon'
    )

    # 6. EXTRACT & SAVE RESULTS
    log_message("Extracting results", "STEP")
    result = adata.uns['rank_genes_groups']
    group = group1
    contrast_label = f"{group1} vs {group2}"
    out_df = pd.DataFrame({
        'gene': result['names'][group],
        'logfc': result['logfoldchanges'][group],
        'pval': result['pvals'][group],
        'pval_adj': result['pvals_adj'][group],
        'score': result['scores'][group],
    })

    # Compute extra columns
    padj_nonzero = out_df['pval_adj'].replace(0, np.nanmin(out_df['pval_adj'][out_df['pval_adj'] > 0]))
    out_df['minlog10padj'] = -np.log10(padj_nonzero)
    out_df['regulationDirection'] = np.select(
        [(out_df['logfc'] >= logfc_cutoff) & (out_df['pval_adj'] <= padj_cutoff), 
            (out_df['logfc'] <= -logfc_cutoff) & (out_df['pval_adj'] <= padj_cutoff)],
        ['Up', 'Down'],
        default='NS'
    )
    out_df.insert(0, 'Contrast', contrast_label)

    # Save all results
    out_df.to_csv(output, index=False)
    log_message(f"All differential expression results saved to {output}", "DONE")

    # Filter by cutoffs
    filtered_df = out_df[
        (out_df['pval_adj'] <= padj_cutoff) &
        (out_df['logfc'].abs() >= logfc_cutoff)
    ]
    filtered_output = output.replace('.csv', '_filtered.csv')
    filtered_df.to_csv(filtered_output, index=False)
    log_message(f"Filtered results (adj_p < {padj_cutoff}, |logfc| > {logfc_cutoff}) saved to {filtered_output}", "DONE")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Differential expression using Scanpy (Wilcoxon, single-cell)")
    parser.add_argument('--expr', required=True, help='Path to expression CSV')
    parser.add_argument('--meta', required=True, help='Path to metadata CSV')
    parser.add_argument('--group1', required=True, help='Numerator group label (e.g., treated)')
    parser.add_argument('--group2', required=True, help='Denominator group label (e.g., control)')
    parser.add_argument('--condition_col', required=True, help='Column in metadata with the group labels (e.g., Genotype or Condition)')
    parser.add_argument('--output', default='de_results.csv', help='Path to output CSV for all genes')
    parser.add_argument('--padj_cutoff', type=float, default=0.05, help='Adjusted p-value threshold')
    parser.add_argument('--logfc_cutoff', type=float, default=1.0, help='Minimum absolute log2 fold change')
    parser.add_argument('--error_output', default='de_errors.csv', help='Path to output CSV for error logs')
    args = parser.parse_args()

    main(
        args.expr,
        args.meta,
        args.group1,
        args.group2,
        args.output,
        args.condition_col,
        args.padj_cutoff,
        args.logfc_cutoff,
        args.error_output
    )
