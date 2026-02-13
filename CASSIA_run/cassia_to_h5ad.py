#!/usr/bin/env python3
"""
Map CASSIA Output to H5AD

Simple script to map CASSIA cluster-level annotations back to individual cells in h5ad.

Usage:
    python map_cassia_to_h5ad.py
"""

import scanpy as sc
import pandas as pd


def map_cassia_to_h5ad(adata, cassia_csv_path, cluster_column='leiden', 
                       output_column='cell_type'):
    """
    Map CASSIA cluster-level annotations back to individual cells in h5ad.
    
    Parameters:
    -----------
    adata : AnnData
        The AnnData object with cluster assignments
    cassia_csv_path : str
        Path to CASSIA output CSV file
    cluster_column : str
        Column name in adata.obs containing cluster labels (default: 'leiden')
    output_column : str
        Column name to create in adata.obs for cell type annotations
        
    Returns:
    --------
    AnnData
        Updated AnnData object with cell type annotations
    """
    print(f"\nMapping CASSIA annotations to individual cells...")
    
    # Load CASSIA results
    cassia_results = pd.read_csv(cassia_csv_path)
    
    print(f"  ✓ Loaded CASSIA results from: {cassia_csv_path}")
    print(f"  ✓ Found annotations for {len(cassia_results)} clusters")
    
    # Display CASSIA output structure
    print(f"\n  CASSIA output columns: {cassia_results.columns.tolist()}")
    
    # Try to identify the cluster column in CASSIA results
    cluster_col_cassia = None
    for col in ['cluster', 'Cluster', 'cluster_id', 'group']:
        if col in cassia_results.columns:
            cluster_col_cassia = col
            break

    # Fallback: if no obvious cluster column, use the first column as cluster IDs
    if cluster_col_cassia is None:
        fallback_col = cassia_results.columns[0]
        print(f"\n  ⚠️  No explicit cluster column found; using first column '{fallback_col}' as cluster ID")
        cluster_col_cassia = fallback_col
    
    # Try to identify cell type columns
    main_type_col = None
    sub_type_col = None
    
    for col in [
        'Predicted Main Cell Type',
        'main_cell_type',
        'cell_type',
        'annotation',
        'Main_Cell_Type'
    ]:
        if col in cassia_results.columns:
            main_type_col = col
            break
    
    for col in [
        'Predicted Sub Cell Types',
        'sub_cell_type',
        'subtype',
        'Sub_Cell_Type',
        'detailed_annotation'
    ]:
        if col in cassia_results.columns:
            sub_type_col = col
            break
    
    if main_type_col is None:
        raise ValueError(f"Cannot find cell type column in CASSIA results. "
                        f"Available columns: {cassia_results.columns.tolist()}")
    
    print(f"\n  Using CASSIA columns:")
    print(f"    - Cluster: '{cluster_col_cassia}'")
    print(f"    - Main cell type: '{main_type_col}'")
    if sub_type_col:
        print(f"    - Sub cell type: '{sub_type_col}'")
    
    # Create mapping dictionary from cluster to cell type
    cluster_to_celltype = {}
    cluster_to_subtype = {}
    
    for _, row in cassia_results.iterrows():
        cluster_id = str(row[cluster_col_cassia])
        
        if pd.notna(row[main_type_col]):
            cluster_to_celltype[cluster_id] = row[main_type_col]
        
        if sub_type_col and pd.notna(row[sub_type_col]):
            cluster_to_subtype[cluster_id] = row[sub_type_col]
    
    # Map to individual cells
    print(f"\n  Mapping clusters to cells using '{cluster_column}' column...")
    adata.obs[output_column] = adata.obs[cluster_column].astype(str).map(cluster_to_celltype)
    
    if cluster_to_subtype:
        adata.obs[f"{output_column}_subtype"] = adata.obs[cluster_column].astype(str).map(cluster_to_subtype)
        print(f"  ✓ Created '{output_column}_subtype' column")
    
    # Check for unmapped cells
    unmapped = adata.obs[output_column].isna().sum()
    if unmapped > 0:
        print(f"\n  ⚠️  Warning: {unmapped} cells could not be mapped")
        print(f"      This might be due to cluster ID mismatch between h5ad and CASSIA CSV")
    else:
        print(f"\n  ✓ Successfully mapped all {adata.n_obs} cells")
    
    # Summary statistics
    print(f"\n  Annotation summary:")
    print(f"    - Total cells: {adata.n_obs}")
    print(f"    - Unique cell types: {adata.obs[output_column].nunique()}")
    if cluster_to_subtype:
        print(f"    - Unique subtypes: {adata.obs[f'{output_column}_subtype'].nunique()}")
    
    print(f"\n  Cell type distribution:")
    cell_type_counts = adata.obs[output_column].value_counts()
    for cell_type, count in cell_type_counts.items():
        print(f"    {cell_type}: {count}")
    
    return adata


def main():
    """
    Main function to map CASSIA output to h5ad.
    """
    # ============================================================================
    # CONFIGURATION - EDIT THESE SETTINGS
    # ============================================================================
    
    # Input files
    INPUT_H5AD = "your_clustered_data.h5ad"  # Your h5ad with clusters
    CASSIA_CSV = "cassia_results_scored.csv"  # CASSIA output CSV
    
    # Column names
    CLUSTER_COLUMN = "leiden"  # Column in h5ad with cluster labels
    OUTPUT_COLUMN = "cell_type"  # Column name for cell type annotations
    
    # Output file
    OUTPUT_H5AD = "annotated_data.h5ad"
    
    # ============================================================================
    # EXECUTION
    # ============================================================================
    
    print("="*70)
    print("CASSIA to H5AD Mapper")
    print("="*70)
    
    # Load h5ad file
    print(f"\nLoading h5ad from: {INPUT_H5AD}")
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"  ✓ Loaded {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Check if cluster column exists
    if CLUSTER_COLUMN not in adata.obs.columns:
        print(f"\n❌ Error: Cluster column '{CLUSTER_COLUMN}' not found in h5ad")
        print(f"   Available columns: {adata.obs.columns.tolist()}")
        return
    
    print(f"  ✓ Found cluster column: '{CLUSTER_COLUMN}'")
    print(f"  ✓ Number of clusters: {adata.obs[CLUSTER_COLUMN].nunique()}")
    
    # Map CASSIA annotations
    adata = map_cassia_to_h5ad(
        adata=adata,
        cassia_csv_path=CASSIA_CSV,
        cluster_column=CLUSTER_COLUMN,
        output_column=OUTPUT_COLUMN
    )
    
    # Save annotated h5ad
    print(f"\nSaving annotated h5ad to: {OUTPUT_H5AD}")
    adata.write(OUTPUT_H5AD)
    
    print("\n" + "="*70)
    print("Mapping completed successfully!")
    print("="*70)
    print(f"\nOutput file: {OUTPUT_H5AD}")
    print(f"Cell type column: '{OUTPUT_COLUMN}'")
    
    if f"{OUTPUT_COLUMN}_subtype" in adata.obs.columns:
        print(f"Subtype column: '{OUTPUT_COLUMN}_subtype'")
    
    print("\n✓ Ready for evaluation!")
    print()


if __name__ == "__main__":
    main()