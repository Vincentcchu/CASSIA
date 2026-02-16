#!/usr/bin/env python3
"""
CASSIA Analysis Pipeline for H5AD Data with Performance Tracking

This script processes single-cell RNA-seq data in h5ad format and runs
CASSIA automated cell type annotation with comprehensive tracking of
time, token usage, and costs.

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
import time


def _load_api_key(config_path: str = "api_keys.json", provider: str = "openai") -> str:
    """Load provider-specific API key from a simple JSON file."""
    with open(config_path, "r", encoding="utf-8") as handle:
        data = json.load(handle)

    key_map = {
        "openai": "openai",
        "openrouter": "openrouter",
        "anthropic": "anthropic",
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
    tuple: (AnnData, preprocessing_time)
        Processed AnnData object with clusters and time taken
    """
    print("Starting preprocessing...")
    start_time = time.time()
    
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
    # Copy the Leiden labels to a CASSIA-friendly column name
    adata.obs['cluster'] = adata.obs['leiden']
    
    # Compute UMAP for visualization
    print("  - Computing UMAP...")
    sc.tl.umap(adata)
    
    preprocessing_time = time.time() - start_time
    
    print(f"✓ Preprocessing complete! Found {len(adata.obs['leiden'].unique())} clusters")
    print(f"  Time taken: {preprocessing_time:.2f} seconds ({preprocessing_time/60:.2f} minutes)")
    
    return adata, preprocessing_time


def find_marker_genes(adata, groupby='cluster', method='wilcoxon'):
    """
    Find marker genes for each cluster and format for CASSIA.
    
    CASSIA expects these columns:
    - cluster: Cluster identifier for the gene
    - gene: Gene symbol
    - avg_log2FC: Average log2 fold change
    - pct.1: Percentage of cells in cluster expressing the gene
    - pct.2: Percentage of cells outside cluster expressing the gene
    - p_val, p_val_adj: (Optional) p-values and adjusted p-values
    
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
        Marker genes dataframe formatted for CASSIA
    """
    print("\nFinding marker genes...")
    
    # Compute marker genes with percentage information (pts=True)
    sc.tl.rank_genes_groups(adata, groupby, method=method, pts=True)
    
    # Extract as dataframe
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    
    print(f"  Raw scanpy columns: {list(markers.columns)}")
    
    # Rename columns to match CASSIA's expected format
    # scanpy outputs: 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj', 
    #                 'group', 'pct_nz_group', 'pct_nz_reference'
    column_mapping = {
        'group': 'cluster',              # Cluster identifier
        'names': 'gene',                 # Gene symbol
        'logfoldchanges': 'avg_log2FC',  # Log2 fold change
        'pct_nz_group': 'pct.1',         # Percentage in cluster
        'pct_nz_reference': 'pct.2',     # Percentage outside cluster
        'pvals': 'p_val',                # P-values
        'pvals_adj': 'p_val_adj'         # Adjusted p-values
    }
    
    # Rename columns
    markers = markers.rename(columns=column_mapping)
    
    # Select only the columns CASSIA expects (in the correct order)
    cassia_columns = ['cluster', 'gene', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val', 'p_val_adj']
    
    # Keep only columns that exist
    available_columns = [col for col in cassia_columns if col in markers.columns]
    markers = markers[available_columns]
    
    print(f"\n✓ Formatted marker genes for CASSIA")
    print(f"  Clusters: {len(markers['cluster'].unique())}")
    print(f"  Total marker records: {len(markers)}")
    print(f"  Columns: {list(markers.columns)}")
    print(f"\n  Sample data:")
    print(markers.head(10).to_string())
    
    return markers


def main():
    """
    Main function to run the complete pipeline with performance tracking.
    """
    # ============================================================================
    # CONFIGURATION - EDIT THESE SETTINGS
    # ============================================================================
    
    # Input file
    INPUT_H5AD = "/cs/student/projects2/aisd/2024/shekchu/projects/data/brain/h5ad_unlabelled/Data_Choudhury2022_Brain.h5ad"  # Change this to your h5ad file path
    
    # CASSIA settings
    PROVIDER = "openai"        # "openrouter", "openai", or "anthropic"
    API_KEY = _load_api_key(provider=PROVIDER)  # Reads from api_keys.json
    MODEL = "gpt-5.1"  # Model to use
    
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
    
    # Cost tracking settings (GPT-5.1 pricing)
    INPUT_PRICE_PER_MILLION = 1.25  # $1.25 per 1M input tokens
    OUTPUT_PRICE_PER_MILLION = 10.0  # $10 per 1M output tokens
    
    # Optional: save intermediate results
    SAVE_MARKERS = True            # Save marker genes to CSV
    SAVE_PROCESSED_H5AD = True     # Save processed h5ad file
    
    # ============================================================================
    # PIPELINE EXECUTION
    # ============================================================================
    
    # Import token tracker
    from CASSIA.token_tracker import get_tracker, reset_tracker
    
    # Reset and start tracking
    reset_tracker()
    tracker = get_tracker()
    tracker.start_timer()
    
    print("="*70)
    print("CASSIA Analysis Pipeline with Performance Tracking")
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
    
    # Preprocess and cluster (with timing)
    adata, preprocessing_time = preprocess_and_cluster(
        adata, 
        resolution=CLUSTERING_RESOLUTION,
        n_pcs=N_PCS,
        n_neighbors=N_NEIGHBORS
    )
    
    # Record preprocessing time
    tracker.record_preprocessing_time(preprocessing_time)
    
    # Save processed h5ad if requested
    if SAVE_PROCESSED_H5AD:
        processed_file = f"{OUTPUT_NAME}_processed.h5ad"
        print(f"\nSaving processed data to: {processed_file}")
        adata.write(processed_file)
    
    # Find marker genes
    markers = find_marker_genes(adata, groupby='cluster', method='wilcoxon')
    
    # Save markers if requested
    if SAVE_MARKERS:
        OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        print(f"\nOutput directory: {os.path.abspath(OUTPUT_DIR)}")
        markers_file = os.path.join(OUTPUT_DIR, f"{OUTPUT_NAME}_markers.csv")
        print(f"\nSaving marker genes to: {markers_file}")
        markers.to_csv(markers_file, index=False)
        
        # Use the saved markers file for CASSIA
        marker_path = markers_file
    else:
        # If not saving, create temporary markers file
        marker_path = f"{OUTPUT_NAME}_markers_temp.csv"
        markers.to_csv(marker_path, index=False)
    
    # ============================================================================
    # RUN CASSIA PIPELINE
    # ============================================================================
    
    print("\n" + "="*70)
    print("Running CASSIA Cell Type Annotation")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"  Provider: {PROVIDER}")
    print(f"  Model: {MODEL}")
    print(f"  Tissue: {TISSUE}")
    print(f"  Species: {SPECIES}")
    print(f"  Clusters: {len(markers['cluster'].unique())}")
    print(f"  Max workers: {MAX_WORKERS}")
    
    # Set API key for CASSIA
    print(f"\nSetting API key for {PROVIDER}...")
    CASSIA.set_api_key(API_KEY, provider=PROVIDER)
    
    # Record annotation start time
    tracker.record_annotation_start()
    
    try:
        CASSIA.runCASSIA_pipeline(
            output_file_name=OUTPUT_NAME,
            tissue=TISSUE,
            species=SPECIES,
            marker=marker_path,
            max_workers=MAX_WORKERS,
            annotation_model=MODEL,
            annotation_provider=PROVIDER,
            score_model=MODEL,
            score_provider=PROVIDER,
            score_threshold=75,
            annotationboost_model=MODEL,
            annotationboost_provider=PROVIDER
        )
        
        # Record annotation end time
        tracker.record_annotation_end()
        
        print("\n" + "="*70)
        print("✓ CASSIA Analysis Complete!")
        print("="*70)
        print(f"\nResults saved with prefix: {OUTPUT_NAME}")
        
    except Exception as e:
        print(f"\n❌ Error running CASSIA pipeline: {e}")
        raise
    finally:
        # Stop overall timer
        tracker.stop_timer()
        
        # Print performance summary
        tracker.print_summary(
            input_price_per_million=INPUT_PRICE_PER_MILLION,
            output_price_per_million=OUTPUT_PRICE_PER_MILLION
        )
        
        # Save performance summary to file
        summary_file = f"{OUTPUT_NAME}_performance.json"
        tracker.save_summary(
            summary_file,
            input_price_per_million=INPUT_PRICE_PER_MILLION,
            output_price_per_million=OUTPUT_PRICE_PER_MILLION
        )
    
    # Clean up temporary markers file if it was created
    if not SAVE_MARKERS and os.path.exists(marker_path):
        os.remove(marker_path)
        print(f"\nCleaned up temporary file: {marker_path}")


if __name__ == "__main__":
    main()