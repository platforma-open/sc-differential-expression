import argparse
import gc
import os
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

# Fixed column types for the result tables. Parquet carries its own types, so an
# empty table must still declare them, or ptabler cannot concatenate it with a
# non-empty one.
RESULT_SCHEMA = {
    "Contrast": pl.String,
    "gene": pl.String,
    "logfc": pl.Float64,
    "pval": pl.Float64,
    "pval_adj": pl.Float64,
    "score": pl.Float64,
    "minlog10padj": pl.Float64,
    "regulationDirection": pl.String,
}
ERROR_SCHEMA = {"Error": pl.String, "value": pl.String}

def write_table(df, path, schema):
    """Write a pandas frame as Parquet with the given column types.

    NaN becomes null, as an empty CSV field did before, so the p-columns keep
    reading a missing value as NA.
    """
    table = pl.DataFrame({c: df[c].tolist() for c in schema}, schema=schema, strict=False)
    float_columns = [c for c, t in schema.items() if t == pl.Float64]
    table = table.with_columns([pl.col(c).fill_nan(None) for c in float_columns])
    table.write_parquet(path)

def filtered_path(output):
    root, ext = os.path.splitext(output)
    return f"{root}_filtered{ext}"

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
        df_empty = pd.DataFrame(columns=list(RESULT_SCHEMA))
        write_table(df_empty, output, RESULT_SCHEMA)
        write_table(df_empty, filtered_path(output), RESULT_SCHEMA)

        # output error logs table
        df_error = pd.DataFrame({"Error": errorLogs, "value": errorLogs})
        write_table(df_error, error_output, ERROR_SCHEMA)
        return
    else:
        df_error = pd.DataFrame(columns=list(ERROR_SCHEMA))
        write_table(df_error, error_output, ERROR_SCHEMA)

    # 2. LOAD EXPRESSION DATA WITH POLARS & CATEGORICAL OPTIMIZATION
    log_message("Loading expression data with Polars optimization", "STEP")

    # Peek schema for flexible headers and identify columns for Categorical casting
    temp_scan = pl.scan_parquet(expr_path)
    column_names = set(temp_scan.collect_schema().names())

    # The repeated string columns are held as categoricals. Parquet carries its own
    # dtypes, so the cast rides in the scan plan instead of a read-time override.
    categorical_candidates = ["Sample", "Ensembl Id", "Cell Barcode", "Cell ID"]
    categorical_columns = [c for c in categorical_candidates if c in column_names]
    df_expr_pl = temp_scan.with_columns(
        [pl.col(c).cast(pl.Categorical) for c in categorical_columns]
    ).collect()
    log_message(f"Loaded {df_expr_pl.height} rows x {df_expr_pl.width} columns", "DONE")

    # Normalize header labels to support migration from 'Cell Barcode' to 'Cell ID'
    if 'Cell Barcode' not in df_expr_pl.columns and 'Cell ID' in df_expr_pl.columns:
        df_expr_pl = df_expr_pl.rename({'Cell ID': 'Cell Barcode'})

    SEPARATOR = '|||'

    # 3. CONSTRUCT SPARSE MATRIX FROM CATEGORICAL CODES
    log_message("Creating sparse matrix from long format data", "STEP")

    # Everything below works on the categoricals' integer codes, not on their
    # strings. Building a per-row "Sample|||Cell Barcode" string column turns
    # both categoricals back into full strings for every row, which runs out of
    # memory on large count matrices. Codes take 4 bytes per row, and the string
    # work happens once per distinct value instead of once per row.
    sample_codes = np.array(df_expr_pl['Sample'].to_physical().to_numpy(), copy=True)
    barcode_codes = np.array(df_expr_pl['Cell Barcode'].to_physical().to_numpy(), copy=True)
    gene_codes = np.array(df_expr_pl['Ensembl Id'].to_physical().to_numpy(), copy=True)
    # Use float32 for expression values (Scanpy standard)
    expression_values = np.array(
        df_expr_pl['Raw gene expression'].cast(pl.Float32).to_numpy(), copy=True)

    sample_categories = df_expr_pl['Sample'].cat.get_categories().to_numpy()
    barcode_categories = df_expr_pl['Cell Barcode'].cat.get_categories().to_numpy()
    gene_categories = df_expr_pl['Ensembl Id'].cat.get_categories().to_numpy()

    # Release the long frame before the remapping stage; it is the largest allocation.
    del df_expr_pl
    gc.collect()

    # Gene axis: sorted unique Ensembl Ids, and each row's index into them.
    used_gene_phys, gene_idx = np.unique(gene_codes, return_inverse=True)
    del gene_codes
    unique_gene_ids, gene_ranks = np.unique(gene_categories[used_gene_phys], return_inverse=True)
    col_codes = gene_ranks[gene_idx].astype(np.int32)
    del gene_idx, gene_ranks, used_gene_phys
    gc.collect()

    # Cell axis: pair the two categorical codes into one integer key, then build
    # the label strings only for the distinct pairs that actually occur.
    n_barcodes = len(barcode_categories)
    pair_codes = sample_codes.astype(np.int64) * n_barcodes + barcode_codes.astype(np.int64)
    del sample_codes, barcode_codes
    gc.collect()

    used_pairs, pair_idx = np.unique(pair_codes, return_inverse=True)
    del pair_codes
    gc.collect()

    pair_labels = np.array(
        [f"{sample_categories[p // n_barcodes]}{SEPARATOR}{barcode_categories[p % n_barcodes]}"
         for p in used_pairs],
        dtype=object)
    unique_cell_ids, cell_ranks = np.unique(pair_labels, return_inverse=True)
    row_codes = cell_ranks[pair_idx].astype(np.int32)
    del pair_idx, cell_ranks, used_pairs, pair_labels
    gc.collect()

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
    write_table(out_df, output, RESULT_SCHEMA)
    log_message(f"All differential expression results saved to {output}", "DONE")

    # Filter by cutoffs
    filtered_df = out_df[
        (out_df['pval_adj'] <= padj_cutoff) &
        (out_df['logfc'].abs() >= logfc_cutoff)
    ]
    filtered_output = filtered_path(output)
    write_table(filtered_df, filtered_output, RESULT_SCHEMA)
    log_message(f"Filtered results (adj_p < {padj_cutoff}, |logfc| > {logfc_cutoff}) saved to {filtered_output}", "DONE")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Differential expression using Scanpy (Wilcoxon, single-cell)")
    parser.add_argument('--expr', required=True, help='Path to long-format raw counts Parquet')
    parser.add_argument('--meta', required=True, help='Path to metadata CSV')
    parser.add_argument('--group1', required=True, help='Numerator group label (e.g., treated)')
    parser.add_argument('--group2', required=True, help='Denominator group label (e.g., control)')
    parser.add_argument('--condition_col', required=True, help='Column in metadata with the group labels (e.g., Genotype or Condition)')
    parser.add_argument('--output', default='de_results.parquet', help='Path to output Parquet for all genes')
    parser.add_argument('--padj_cutoff', type=float, default=0.05, help='Adjusted p-value threshold')
    parser.add_argument('--logfc_cutoff', type=float, default=1.0, help='Minimum absolute log2 fold change')
    parser.add_argument('--error_output', default='de_errors.parquet', help='Path to output Parquet for error logs')
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
