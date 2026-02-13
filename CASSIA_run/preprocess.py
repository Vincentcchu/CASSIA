#!/usr/bin/env python3
"""
CASSIA Analysis Pipeline for H5AD Data

This script processes single-cell RNA-seq data in h5ad format and runs
CASSIA automated cell type annotation.

Usage:
    python cassia_analysis.py

Requirements:
    pip install scanpy CASSIA pandas numpy
"""

import scanpy as sc
import pandas as pd
import CASSIA
import argparse
import os
import json


def _load_api_key(config_path: str = "../api_keys.json", provider: str = "openai") -> str:
    """Load provider-specific API key from a simple JSON file."""
    with open(config_path, "r", encoding="utf-8") as handle:
        data = json.load(handle)

    key_map = {
        "openai": "openai_api_key",
        "openrouter": "openrouter_api_key",
        "anthropic": "anthropic_api_key",
    }
    key_field = key_map.get(provider.lower())
    if not key_field:
        raise ValueError(f"Unknown provider '{provider}'. Expected one of: {', '.join(key_map)}")
    if key_field not in data:
        raise KeyError(f"Key '{key_field}' missing in {config_path}")
    return data[key_field]

def preprocess_and_cluster(adata, resolution=0.5, n_pcs=40, n_neighbors=10):
    """
    Perform standard preprocessing and clustering on AnnData object.
    
    Parameters:
    -----------
    adata : AnnData
        The input AnnData object
    resolution : float
        Resolution parameter for Leiden clustering
    n_pcs : int
        Number of principal components to compute
    n_neighbors : int
        Number of neighbors for graph construction
        
    Returns:
    --------
    AnnData
        Processed AnnData object with clusters
    """
    print("Starting preprocessing...")
    
    # Basic QC and filtering
    print("  - Filtering cells and genes...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Normalize and log transform
    print("  - Normalizing data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Store raw data
    adata.raw = adata
    
    # Find highly variable genes
    print("  - Identifying highly variable genes...")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    
    # Scale data
    print("  - Scaling data...")
    sc.pp.scale(adata, max_value=10)
    
    # Run PCA
    print("  - Computing PCA...")
    sc.tl.pca(adata, svd_solver='arpack')
    
    # Compute neighborhood graph
    print("  - Building neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # Cluster the cells
    print(f"  - Clustering with Leiden (resolution={resolution})...")
    sc.tl.leiden(adata, resolution=resolution)
    
    # Compute UMAP for visualization
    print("  - Computing UMAP...")
    sc.tl.umap(adata)
    
    print(f"✓ Preprocessing complete! Found {len(adata.obs['leiden'].unique())} clusters")
    
    return adata


def find_marker_genes(adata, groupby='leiden', method='wilcoxon'):
    """
    Find marker genes for each cluster.
    
    Parameters:
    -----------
    adata : AnnData
        Processed AnnData object with clusters
    groupby : str
        Column name in adata.obs containing cluster labels
    method : str
        Method for differential expression ('wilcoxon', 't-test', 'logreg')
        
    Returns:
    --------
    DataFrame
        Marker genes dataframe for CASSIA
    """
    print("\nFinding marker genes...")
    
    # Compute marker genes
    sc.tl.rank_genes_groups(adata, groupby, method=method)
    
    # Extract as dataframe
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    
    print(f"✓ Found marker genes for {len(markers['group'].unique())} clusters")
    print(f"  Total marker records: {len(markers)}")
    
    return markers


def main():
    """
    Main function to run the complete pipeline.
    """
    # ============================================================================
    # CONFIGURATION - EDIT THESE SETTINGS
    # ============================================================================
    
    # Input file
    INPUT_H5AD = "/cs/student/projects2/aisd/2024/shekchu/projects/data/brain/h5ad_unlabelled/Data_Choudhury2022_Brain.h5ad"  # Change this to your h5ad file path
    
    # CASSIA settings
    PROVIDER = "openai"        # "openrouter", "openai", or "anthropic"
    API_KEY = _load_api_key(provider=PROVIDER)  # Reads from api_keys.json
    MODEL = "gpt-4o-2024-08-06"  # Model to use
    
    # Biological context
    TISSUE = "Brain"               # Your tissue type
    SPECIES = "Human"              # "Human" or "Mouse"
    
    # Output settings
    OUTPUT_NAME = "cassia_results"  # Prefix for output files
    
    # Analysis parameters
    CLUSTERING_RESOLUTION = 0.5    # Leiden clustering resolution (higher = more clusters)
    N_PCS = 40                     # Number of principal components
    N_NEIGHBORS = 10               # Number of neighbors for graph
    MAX_WORKERS = 4                # Parallel workers for CASSIA
    
    # Optional: save intermediate results
    SAVE_MARKERS = True            # Save marker genes to CSV
    SAVE_PROCESSED_H5AD = True     # Save processed h5ad file
    
    # ============================================================================
    # PIPELINE EXECUTION
    # ============================================================================
    
    print("="*70)
    print("CASSIA Analysis Pipeline")
    print("="*70)
    
    # Check if input file exists
    if not os.path.exists(INPUT_H5AD):
        print(f"\n❌ Error: Input file '{INPUT_H5AD}' not found!")
        print("Please update the INPUT_H5AD variable with your h5ad file path.")
        return
    
    # Load h5ad file
    print(f"\nLoading data from: {INPUT_H5AD}")
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"✓ Loaded {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Preprocess and cluster
    adata = preprocess_and_cluster(
        adata, 
        resolution=CLUSTERING_RESOLUTION,
        n_pcs=N_PCS,
        n_neighbors=N_NEIGHBORS
    )
    
    # Save processed h5ad if requested
    if SAVE_PROCESSED_H5AD:
        processed_file = f"{OUTPUT_NAME}_processed.h5ad"
        print(f"\nSaving processed data to: {processed_file}")
        adata.write(processed_file)
    
    # Find marker genes
    markers = find_marker_genes(adata, groupby='leiden', method='wilcoxon')
    
    # Save markers if requested
    if SAVE_MARKERS:
        OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        print(f"\nOutput directory: {os.path.abspath(OUTPUT_DIR)}")
        markers_file = os.path.join(OUTPUT_DIR, f"{OUTPUT_NAME}_markers.csv")
        print(f"\nSaving marker genes to: {markers_file}")
        markers.to_csv(markers_file, index=False)


if __name__ == "__main__":
    main()