#!/usr/bin/env python3
"""
Batch processing script for running CASSIA on multiple h5ad files.

This script can process:
- All tissues automatically
- A specific tissue type
- Specific files within a tissue

Output Structure:
    batch_results/run_{tissue}_{YYYYMMDD_HHMMSS}/  # Tissue-labeled timestamped run directory
    ├── batch_summary.txt                          # Human-readable summary
    ├── batch_summary.csv                          # Detailed results table
    ├── token_tracking.json                        # Token usage details
    └── {tissue_name}/                             # Per-tissue folders
        └── {dataset_name}/                        # Per-dataset folders
            ├── {dataset_name}_markers.csv         # Marker genes
            ├── {dataset_name}_annotated.h5ad      # Annotated data
            ├── {dataset_name}_annotation.txt      # CASSIA text output
            └── {dataset_name}_summary.txt         # Individual file summary

Naming Convention:
    - Run directories include tissue name(s):
      * Single tissue: run_brain_20260216_151425
      * Multiple tissues: run_brain_breast_20260216_151425
      * All tissues: run_all_20260216_151425
    - Dataset name is extracted from the h5ad filename (e.g., "Data_Bassez2021_Breast")
    - All output files for a dataset use this consistent prefix
    - Timestamped to avoid conflicts between runs

Usage examples:
    # Process all tissues
    python batch_process_all_tissues.py
    
    # Process only brain tissue
    python batch_process_all_tissues.py --tissue brain
    
    # Process specific files in breast tissue
    python batch_process_all_tissues.py --tissue breast --files "Data_Bassez*.h5ad"
    
    # Process multiple specific tissues
    python batch_process_all_tissues.py --tissues brain breast colorectal
    
    # Preview without processing
    python batch_process_all_tissues.py --preview
"""

import os
import sys
import json
import time
import glob
import pandas as pd
import scanpy as sc
import argparse
from pathlib import Path
from datetime import datetime
import CASSIA

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(__file__))

# Import functions from the main workflow
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
    """Perform standard preprocessing and clustering on AnnData object."""
    print("  Starting preprocessing...")
    
    # Basic QC and filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Normalize and log transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Store raw data
    adata.raw = adata
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    
    # Scale data
    sc.pp.scale(adata, max_value=10)
    
    # Run PCA
    sc.tl.pca(adata, svd_solver='arpack')
    
    # Compute neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # Cluster the cells
    sc.tl.leiden(adata, resolution=resolution)
    adata.obs['cluster'] = adata.obs['leiden']
    
    # Compute UMAP
    sc.tl.umap(adata)
    
    print(f"  ✓ Preprocessing complete! Found {len(adata.obs['leiden'].unique())} clusters")
    
    return adata


def find_marker_genes(adata, groupby='cluster', method='wilcoxon'):
    """Find marker genes for each cluster and format for CASSIA."""
    print("  Finding marker genes...")
    
    # Compute marker genes with percentage information
    sc.tl.rank_genes_groups(adata, groupby, method=method, pts=True)
    
    # Extract as dataframe
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    
    # Rename columns to match CASSIA's expected format
    column_mapping = {
        'group': 'cluster',
        'names': 'gene',
        'logfoldchanges': 'avg_log2FC',
        'pct_nz_group': 'pct.1',
        'pct_nz_reference': 'pct.2',
        'pvals': 'p_val',
        'pvals_adj': 'p_val_adj'
    }
    
    markers = markers.rename(columns=column_mapping)
    
    # Select only the columns CASSIA expects
    cassia_columns = ['cluster', 'gene', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val', 'p_val_adj']
    available_columns = [col for col in cassia_columns if col in markers.columns]
    markers = markers[available_columns]
    
    print(f"  ✓ Found {len(markers['cluster'].unique())} clusters with markers")
    
    return markers


def discover_h5ad_files(data_root: str, tissue_types: list = None):
    """
    Discover all h5ad files organized by tissue type.
    
    Args:
        data_root: Root directory containing tissue folders
        tissue_types: List of specific tissue types to process (None = all)
    
    Returns:
        Dictionary mapping tissue type to list of h5ad file paths
    """
    data_path = Path(data_root)
    files_by_tissue = {}
    
    # If specific tissues provided, use those; otherwise scan all directories
    if tissue_types:
        tissue_dirs = [data_path / tissue for tissue in tissue_types if (data_path / tissue).exists()]
    else:
        tissue_dirs = [d for d in data_path.iterdir() if d.is_dir()]
    
    for tissue_dir in tissue_dirs:
        tissue_name = tissue_dir.name
        
        # Look for h5ad files in the h5ad_unlabelled subdirectory
        unlabelled_dir = tissue_dir / "h5ad_unlabelled"
        
        if unlabelled_dir.exists():
            h5ad_files = sorted(unlabelled_dir.glob("*.h5ad"))
            if h5ad_files:
                files_by_tissue[tissue_name] = [str(f) for f in h5ad_files]
                print(f"Found {len(h5ad_files)} h5ad files in {tissue_name}/")
    
    return files_by_tissue


def process_single_file(
    h5ad_path: str,
    tissue_name: str,
    species: str,
    provider: str,
    api_key: str,
    model: str,
    output_base_dir: str,
    config: dict
):
    """
    Process a single h5ad file through the CASSIA pipeline.
    
    Output naming convention:
    - {dataset_name}_markers.csv          : Marker genes for each cluster
    - {dataset_name}_annotated.h5ad       : Annotated AnnData object
    - {dataset_name}_annotation.txt       : CASSIA annotation text output
    - {dataset_name}_summary.txt          : Processing summary and statistics
    
    Returns:
        Dictionary with processing results and statistics
    """
    file_start_time = time.time()
    
    # Extract dataset name from filename (e.g., "Data_Bassez2021_Breast")
    dataset_name = Path(h5ad_path).stem
    
    print("\n" + "="*70)
    print(f"Processing: {dataset_name}")
    print(f"Tissue: {tissue_name}")
    print(f"File: {h5ad_path}")
    print("="*70)
    
    try:
        # Load data
        print("\n[1/4] Loading data...")
        adata = sc.read_h5ad(h5ad_path)
        print(f"✓ Loaded {adata.n_obs} cells × {adata.n_vars} genes")
        
        # Preprocess and cluster
        print("\n[2/4] Preprocessing and clustering...")
        preprocessing_start = time.time()
        adata = preprocess_and_cluster(
            adata,
            resolution=config['clustering_resolution'],
            n_pcs=config['n_pcs'],
            n_neighbors=config['n_neighbors']
        )
        preprocessing_time = time.time() - preprocessing_start
        
        # Find marker genes
        print("\n[3/4] Finding marker genes...")
        markers = find_marker_genes(adata, groupby='cluster', method='wilcoxon')
        
        # Create output directory structure: output_base_dir/tissue_name/dataset_name/
        output_dir = Path(output_base_dir) / tissue_name / dataset_name
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"✓ Output directory: {output_dir}")
        
        # Save markers with clear naming convention
        markers_file = output_dir / f"{dataset_name}_markers.csv"
        markers.to_csv(markers_file, index=False)
        print(f"✓ Markers saved: {markers_file.name}")
        
        # Run CASSIA annotation
        print("\n[4/4] Running CASSIA annotation...")
        output_name = str(output_dir / dataset_name)
        
        annotation_start = time.time()
        CASSIA.runCASSIA_pipeline(
            output_file_name=output_name,
            tissue=tissue_name,
            species=species,
            marker_path=str(markers_file),
            max_workers=config['max_workers'],
            annotation_model=model,
            annotation_provider=provider,
            score_model=model,
            score_provider=provider,
            score_threshold=config['score_threshold'],
            annotationboost_model=model,
            annotationboost_provider=provider
        )
        annotation_time = time.time() - annotation_start
        
        file_total_time = time.time() - file_start_time
        
        print(f"\n✓ {dataset_name} completed successfully!")
        print(f"  Preprocessing: {preprocessing_time:.2f}s")
        print(f"  Annotation: {annotation_time:.2f}s")
        print(f"  Total: {file_total_time:.2f}s")
        
        return {
            'status': 'success',
            'dataset': dataset_name,
            'tissue': tissue_name,
            'file_path': h5ad_path,
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'n_clusters': len(adata.obs['leiden'].unique()),
            'preprocessing_time': preprocessing_time,
            'annotation_time': annotation_time,
            'total_time': file_total_time,
            'output_dir': str(output_dir)
        }
        
    except Exception as e:
        print(f"\n❌ Error processing {dataset_name}: {e}")
        import traceback
        traceback.print_exc()
        
        return {
            'status': 'failed',
            'dataset': dataset_name,
            'tissue': tissue_name,
            'file_path': h5ad_path,
            'error': str(e),
            'total_time': time.time() - file_start_time
        }


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Batch process h5ad files with CASSIA',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all tissues
  %(prog)s
  
  # Process only brain tissue
  %(prog)s --tissue brain
  
  # Process multiple specific tissues
  %(prog)s --tissues brain breast colorectal
  
  # Process specific files in a tissue
  %(prog)s --tissue breast --files "Data_Bassez*.h5ad"
  
  # Preview without processing
  %(prog)s --preview
  
  # Use different model
  %(prog)s --provider openrouter --model "google/gemini-2.5-flash"
        """
    )
    
    # Data selection
    parser.add_argument('--tissue', type=str,
                        help='Process a single tissue type (e.g., brain, breast, colorectal)')
    parser.add_argument('--tissues', nargs='+',
                        help='Process multiple specific tissues (e.g., --tissues brain breast)')
    parser.add_argument('--files', default='*.h5ad',
                        help='File pattern to match (default: *.h5ad)')
    parser.add_argument('--data-root', 
                        default='/cs/student/projects2/aisd/2024/shekchu/projects/data',
                        help='Root data directory (default: %(default)s)')
    
    # Model settings
    parser.add_argument('--provider', default='openai',
                        choices=['openai', 'anthropic', 'openrouter'],
                        help='LLM provider (default: %(default)s)')
    parser.add_argument('--model', default='gpt-4o-2024-08-06',
                        help='Model name (default: %(default)s)')
    parser.add_argument('--species', default='Human',
                        help='Species (default: %(default)s)')
    
    # Processing parameters
    parser.add_argument('--resolution', type=float, default=0.5,
                        help='Clustering resolution (default: %(default)s)')
    parser.add_argument('--n-pcs', type=int, default=40,
                        help='Number of principal components (default: %(default)s)')
    parser.add_argument('--n-neighbors', type=int, default=10,
                        help='Number of neighbors (default: %(default)s)')
    parser.add_argument('--max-workers', type=int, default=4,
                        help='Maximum parallel workers (default: %(default)s)')
    parser.add_argument('--score-threshold', type=float, default=75,
                        help='Score threshold for boost annotation (default: %(default)s)')
    
    # Output
    parser.add_argument('--output-dir', default='./batch_results',
                        help='Output directory (default: %(default)s)')
    
    # Cost tracking
    parser.add_argument('--input-cost', type=float, default=1.25,
                        help='Input cost per 1M tokens in USD (default: %(default)s)')
    parser.add_argument('--output-cost', type=float, default=10.0,
                        help='Output cost per 1M tokens in USD (default: %(default)s)')
    
    # Mode
    parser.add_argument('--preview', action='store_true',
                        help='Preview files without processing')
    
    return parser.parse_args()


def main():
    """Main batch processing function."""
    
    # Parse command-line arguments
    args = parse_arguments()
    
    # ============================================================================
    # CONFIGURATION FROM ARGUMENTS
    # ============================================================================
    
    # Data location
    DATA_ROOT = args.data_root
    
    # Determine which tissues to process
    if args.tissue:
        TISSUE_TYPES = [args.tissue]
    elif args.tissues:
        TISSUE_TYPES = args.tissues
    else:
        TISSUE_TYPES = None  # Process all tissues
    
    # Model settings (API key loaded later if not preview)
    PROVIDER = args.provider
    MODEL = args.model
    SPECIES = args.species
    
    # Create timestamped output directory for this run with tissue name
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    if TISSUE_TYPES is None:
        tissue_label = "all"
    elif len(TISSUE_TYPES) == 1:
        tissue_label = TISSUE_TYPES[0]
    else:
        tissue_label = "_".join(TISSUE_TYPES)
    
    OUTPUT_BASE_DIR = Path(args.output_dir) / f"run_{tissue_label}_{timestamp}"
    OUTPUT_BASE_DIR.mkdir(parents=True, exist_ok=True)
    
    # Processing parameters
    CONFIG = {
        'clustering_resolution': args.resolution,
        'n_pcs': args.n_pcs,
        'n_neighbors': args.n_neighbors,
        'max_workers': args.max_workers,
        'score_threshold': args.score_threshold,
    }
    
    # Cost tracking
    INPUT_COST_PER_MILLION = args.input_cost
    OUTPUT_COST_PER_MILLION = args.output_cost
    
    # ============================================================================
    # BATCH PROCESSING
    # ============================================================================
    
    batch_start_time = time.time()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    print("="*70)
    print("CASSIA BATCH PROCESSING")
    print("="*70)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Data root: {DATA_ROOT}")
    print(f"Output directory: {OUTPUT_BASE_DIR}")
    print(f"Model: {MODEL}")
    print(f"Provider: {PROVIDER}")
    print(f"Species: {SPECIES}")
    print("="*70)
    
    # Discover all h5ad files
    print("\nDiscovering h5ad files...")
    files_by_tissue = discover_h5ad_files(DATA_ROOT, TISSUE_TYPES)
    
    # Apply file pattern filter if specified
    if args.files != '*.h5ad':
        import fnmatch
        filtered_files_by_tissue = {}
        for tissue, files in files_by_tissue.items():
            filtered = [f for f in files if fnmatch.fnmatch(Path(f).name, args.files)]
            if filtered:
                filtered_files_by_tissue[tissue] = filtered
        files_by_tissue = filtered_files_by_tissue
        print(f"Applied filter: {args.files}")
    
    if not files_by_tissue:
        print("❌ No h5ad files found!")
        if args.tissue:
            print(f"   Make sure {DATA_ROOT}/{args.tissue}/h5ad_unlabelled/ exists")
        return
    
    total_files = sum(len(files) for files in files_by_tissue.values())
    print(f"\n✓ Found {total_files} total files across {len(files_by_tissue)} tissues")
    
    # If preview mode, show files and exit
    if args.preview:
        print("\n" + "="*70)
        print("PREVIEW MODE - Files to be processed:")
        print("="*70)
        for tissue_name, h5ad_files in files_by_tissue.items():
            print(f"\n{tissue_name.upper()} ({len(h5ad_files)} files):")
            for i, h5ad_path in enumerate(h5ad_files, 1):
                file_name = Path(h5ad_path).name
                try:
                    adata = sc.read_h5ad(h5ad_path)
                    print(f"  {i:2d}. {file_name}")
                    print(f"      {adata.n_obs:,} cells × {adata.n_vars:,} genes")
                except Exception as e:
                    print(f"  {i:2d}. {file_name} (could not read: {e})")
        
        # Estimate cost and time
        avg_time_per_file = 8  # minutes
        total_minutes = total_files * avg_time_per_file
        avg_tokens_per_file = 150000
        total_tokens = total_files * avg_tokens_per_file
        input_tokens = total_tokens * 0.75
        output_tokens = total_tokens * 0.25
        cost = (input_tokens / 1_000_000 * INPUT_COST_PER_MILLION) + (output_tokens / 1_000_000 * OUTPUT_COST_PER_MILLION)
        
        print(f"\n{'='*70}")
        print("ESTIMATES:")
        print(f"{'='*70}")
        print(f"Time: ~{total_minutes:.0f} minutes ({total_minutes/60:.1f} hours)")
        print(f"Tokens: ~{total_tokens/1_000_000:.1f}M tokens")
        print(f"Cost: ~${cost:.2f}")
        print(f"\nTo process these files, run without --preview flag")
        print(f"{'='*70}")
        return
    
    # Load API key (only needed for actual processing)
    API_KEY = _load_api_key(provider=PROVIDER)
    
    # Set API key
    CASSIA.set_api_key(API_KEY, provider=PROVIDER)
    
    # Reset global tracker
    tracker = CASSIA.get_tracker()
    tracker.reset()
    tracker.start_tracking()
    
    # Process all files
    results = []
    file_count = 0
    
    for tissue_name, h5ad_files in files_by_tissue.items():
        print(f"\n{'='*70}")
        print(f"Processing tissue: {tissue_name.upper()}")
        print(f"Files to process: {len(h5ad_files)}")
        print(f"{'='*70}")
        
        for h5ad_path in h5ad_files:
            file_count += 1
            print(f"\n[File {file_count}/{total_files}]")
            
            result = process_single_file(
                h5ad_path=h5ad_path,
                tissue_name=tissue_name,
                species=SPECIES,
                provider=PROVIDER,
                api_key=API_KEY,
                model=MODEL,
                output_base_dir=OUTPUT_BASE_DIR,
                config=CONFIG
            )
            results.append(result)
    
    # Stop tracking
    tracker.stop_tracking()
    batch_total_time = time.time() - batch_start_time
    
    # ============================================================================
    # GENERATE SUMMARY REPORTS
    # ============================================================================
    
    print("\n" + "="*70)
    print("BATCH PROCESSING COMPLETE")
    print("="*70)
    
    # Calculate statistics
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'failed']
    
    print(f"\n📊 Summary:")
    print(f"  Total files: {len(results)}")
    print(f"  Successful: {len(successful)}")
    print(f"  Failed: {len(failed)}")
    print(f"  Total time: {batch_total_time/60:.2f} minutes")
    
    if successful:
        avg_time = sum(r['total_time'] for r in successful) / len(successful)
        total_cells = sum(r['n_cells'] for r in successful)
        total_clusters = sum(r['n_clusters'] for r in successful)
        print(f"  Average time per file: {avg_time:.2f}s")
        print(f"  Total cells processed: {total_cells:,}")
        print(f"  Total clusters annotated: {total_clusters}")
    
    # Print token usage summary
    tracker.print_summary(
        input_cost_per_million=INPUT_COST_PER_MILLION,
        output_cost_per_million=OUTPUT_COST_PER_MILLION
    )
    
    # Save detailed results to CSV in the run directory
    results_df = pd.DataFrame(results)
    results_csv = OUTPUT_BASE_DIR / "batch_summary.csv"
    results_df.to_csv(results_csv, index=False)
    print(f"\n📄 Results summary saved to: {results_csv}")
    
    # Save token tracking report in the run directory
    tracker_report = OUTPUT_BASE_DIR / "token_tracking.json"
    tracker.save_report(
        str(tracker_report),
        input_cost_per_million=INPUT_COST_PER_MILLION,
        output_cost_per_million=OUTPUT_COST_PER_MILLION
    )
    print(f"📄 Token tracking saved to: {tracker_report}")
    
    # Save a human-readable summary
    summary_file = OUTPUT_BASE_DIR / "batch_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("CASSIA BATCH PROCESSING SUMMARY\n")
        f.write("="*70 + "\n")
        f.write(f"Run timestamp: {timestamp.replace('_', ' ')}\n")
        f.write(f"Output directory: {OUTPUT_BASE_DIR}\n")
        f.write(f"Model: {MODEL} ({PROVIDER})\n")
        f.write(f"Species: {SPECIES}\n\n")
        
        f.write(f"Files processed: {successful}/{total_files}\n")
        f.write(f"Files failed: {failed_count}\n")
        f.write(f"Total time: {total_time:.2f}s ({total_time/60:.1f} min)\n\n")
        
        f.write(f"Total cells processed: {total_cells:,}\n")
        f.write(f"Total clusters annotated: {total_clusters}\n\n")
        
        # Token usage
        summary = tracker.get_summary()
        f.write("TOKEN USAGE:\n")
        f.write(f"  Input tokens: {summary['total_input_tokens']:,}\n")
        f.write(f"  Output tokens: {summary['total_output_tokens']:,}\n")
        f.write(f"  Total tokens: {summary['total_tokens']:,}\n")
        f.write(f"  Total cost: ${summary['total_cost']:.2f}\n\n")
        
        # Per-file details
        if successful > 0:
            f.write("SUCCESSFUL FILES:\n")
            for r in [x for x in results if x['status'] == 'success']:
                f.write(f"  {r['dataset']} ({r['tissue']})\n")
                f.write(f"    - {r['n_cells']:,} cells, {r['n_clusters']} clusters\n")
                f.write(f"    - Time: {r['total_time']:.1f}s\n")
                f.write(f"    - Output: {r['output_dir']}\n")
        
        if failed_count > 0:
            f.write("\nFAILED FILES:\n")
            for r in [x for x in results if x['status'] == 'failed']:
                f.write(f"  {r['dataset']} ({r['tissue']}): {r['error']}\n")
    
    print(f"📄 Summary report saved to: {summary_file}")
    
    # Print failed files if any
    if failed:
        print(f"\n⚠️  Failed files:")
        for r in failed:
            print(f"  - {r['dataset']} ({r['tissue']}): {r['error']}")
    
    print("\n" + "="*70)
    print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*70)


if __name__ == "__main__":
    main()


# Preview all brain files
# python batch_process_all_tissues.py --tissue brain --preview

# Process specific files
# python batch_process_all_tissues.py --tissue breast --files "Data_Bassez*.h5ad"

# Process multiple tissues
# python batch_process_all_tissues.py --tissues brain breast colorectal

# Process all with custom model
# python batch_process_all_tissues.py --provider openrouter --model "google/gemini-2.5-flash"