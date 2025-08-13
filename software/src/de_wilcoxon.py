import argparse
import pandas as pd
import numpy as np
import scanpy as sc
from anndata import AnnData

def main(expr_path, meta_path, group1, group2, output, condition_col, padj_cutoff, logfc_cutoff):
    # Load expression data
    df_expr = pd.read_csv(expr_path)
    df_meta = pd.read_csv(meta_path)

    # Normalize header labels to support migration from 'Cell Barcode' to 'Cell ID'
    df_expr.columns = [c.strip() for c in df_expr.columns]
    if 'Cell Barcode' not in df_expr.columns and 'Cell ID' in df_expr.columns:
        df_expr = df_expr.rename(columns={'Cell ID': 'Cell Barcode'})

    # Create unique cell IDs
    df_expr['CellID'] = df_expr['Sample'].astype(str) + '_' + df_expr['Cell Barcode']

    # Pivot to cell-by-gene matrix
    df_matrix = df_expr.pivot_table(
        index='CellID',
        columns='Ensembl Id',
        values='Raw gene expression',
        fill_value=0
    )

    # Map metadata to cells
    cell_to_sample = df_expr.drop_duplicates('CellID')[['CellID', 'Sample']].set_index('CellID')
    cell_metadata = cell_to_sample.join(df_meta.set_index('Sample'), on='Sample')
    cell_metadata = cell_metadata.loc[df_matrix.index]

    # Create AnnData object
    adata = AnnData(X=df_matrix.values, obs=cell_metadata, var=pd.DataFrame(index=df_matrix.columns))

    # Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Subset to relevant groups
    adata.obs[condition_col] = adata.obs[condition_col].astype(str)
    adata = adata[adata.obs[condition_col].isin([group1, group2])]

    # Run DE test
    sc.tl.rank_genes_groups(
        adata,
        groupby=condition_col,
        groups=[group1],
        reference=group2,
        method='wilcoxon'
    )

    # Extract results
    result = adata.uns['rank_genes_groups']
    group = group1
    contrast_label = f"{group1}_vs_{group2}"
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
        [out_df['logfc'] > 0, out_df['logfc'] < 0],
        ['Up', 'Down'],
        default='None'
    )
    out_df.insert(0, 'Contrast', contrast_label)

    # Save all results
    out_df.to_csv(output, index=False)
    print(f"All differential expression results saved to {output}")

    # Filter by cutoffs
    filtered_df = out_df[
        (out_df['pval_adj'] < padj_cutoff) &
        (out_df['logfc'].abs() > logfc_cutoff)
    ]
    filtered_output = output.replace('.csv', '_filtered.csv')
    filtered_df.to_csv(filtered_output, index=False)
    print(f"Filtered results (adj_p < {padj_cutoff}, |logfc| > {logfc_cutoff}) saved to {filtered_output}")

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
    args = parser.parse_args()

    main(
        args.expr,
        args.meta,
        args.group1,
        args.group2,
        args.output,
        args.condition_col,
        args.padj_cutoff,
        args.logfc_cutoff
    )
