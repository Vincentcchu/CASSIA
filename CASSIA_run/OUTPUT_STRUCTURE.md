# CASSIA Batch Processing Output Structure

## Directory Organization

Each batch processing run creates a timestamped directory with tissue name(s) containing all outputs:

```
batch_results/
└── run_brain_20260216_151425/                     # Tissue-labeled timestamped run directory
    │
    ├── batch_summary.txt                          # Human-readable summary
    ├── batch_summary.csv                          # Structured results (CSV)
    ├── token_tracking.json                        # API usage details
    │
    ├── brain/                                     # Tissue folder
    │   ├── Data_Choudhury2022_Brain/              # Dataset folder
    │   │   ├── Data_Choudhury2022_Brain_markers.csv
    │   │   ├── Data_Choudhury2022_Brain_annotated.h5ad
    │   │   ├── Data_Choudhury2022_Brain_annotation.txt
    │   │   └── Data_Choudhury2022_Brain_summary.txt
    │   │
    │   ├── Data_Filbin2018_Brain/
    │   │   └── [same structure as above]
    │   └── ...
    │
    └── breast/
        ├── Data_Bassez2021_Breast/
        │   └── [same structure as above]
        └── ...
```

## File Types

### Run-Level Files

Located in: `batch_results/run_YYYYMMDD_HHMMSS/`

| File | Description | Format | Use Case |
|------|-------------|--------|----------|
| `batch_summary.txt` | Overall summary report | Plain text | Quick review, reporting |
| `batch_summary.csv` | Detailed results table | CSV | Analysis in pandas/Excel |
| `token_tracking.json` | API call details | JSON | Cost analysis, debugging |

**batch_summary.txt** contains:
- Processing statistics (files processed, time elapsed)
- Token usage (input/output/total)
- Cost breakdown
- Per-file success/failure status

**batch_summary.csv** includes columns:
- `dataset`, `tissue`, `file_path`
- `n_cells`, `n_genes`, `n_clusters`
- `preprocessing_time`, `annotation_time`, `total_time`
- `status`, `output_dir`

**token_tracking.json** tracks:
- Total API calls
- Input/output tokens per call
- Timestamps and costs

### Dataset-Level Files

Located in: `batch_results/run_YYYYMMDD_HHMMSS/{tissue}/{dataset}/`

| File | Description | Format | Use Case |
|------|-------------|--------|----------|
| `{dataset}_markers.csv` | Marker genes per cluster | CSV | DEG analysis, validation |
| `{dataset}_annotated.h5ad` | Annotated single-cell data | h5ad | Downstream analysis |
| `{dataset}_annotation.txt` | CASSIA text output | Plain text | Review annotations |
| `{dataset}_summary.txt` | Individual processing log | Plain text | Debugging, QC |

**{dataset}_markers.csv** contains:
- Cluster identifiers
- Gene names
- Statistical scores (p-values, log fold changes)
- Used by CASSIA for cell type annotation

**{dataset}_annotated.h5ad** includes:
- Original expression data
- Clustering results (`.obs['leiden']`)
- Cell type annotations (`.obs['cell_type']`)
- UMAP/PCA embeddings (`.obsm`)

## Naming Convention

### Dataset Names
- Extracted from input h5ad filenames
- Example: `Data_Bassez2021_Breast.h5ad` → `Data_Bassez2021_Breast`
- Consistent across all output files

### Run Directories
- Format: `run_{tissue}_{YYYYMMDD_HHMMSS}`
- Examples:
  - Single tissue: `run_brain_20260216_151425`
  - Multiple tissues: `run_brain_breast_20260216_151425`
  - All tissues: `run_all_20260216_151425`
- Prevents overwriting previous runs
- Easy to identify what was processed and when

## Usage Examples

### Python Analysis

```python
import scanpy as sc
import pandas as pd
from pathlib import Path

# Set your run directory (includes tissue name)
run_dir = Path("batch_results/run_brain_20260216_151425")

# Load the summary
summary = pd.read_csv(run_dir / "batch_summary.csv")
print(f"Processed {len(summary)} datasets")

# Load a specific annotated dataset
adata = sc.read_h5ad(
    run_dir / "brain" / "Data_Filbin2018_Brain" / "Data_Filbin2018_Brain_annotated.h5ad"
)
print(adata.obs['cell_type'].value_counts())

# Load marker genes
markers = pd.read_csv(
    run_dir / "brain" / "Data_Filbin2018_Brain" / "Data_Filbin2018_Brain_markers.csv"
)
print(f"Found {len(markers)} marker genes across clusters")
```

### Parse Token Usage

```python
import json

# Load token tracking
with open(run_dir / "token_tracking.json", 'r') as f:
    tracking = json.load(f)

print(f"Total API calls: {tracking['total_calls']}")
print(f"Total tokens: {tracking['total_tokens']:,}")
print(f"Total cost: ${tracking['total_cost']:.2f}")

# Per-call breakdown
for call in tracking['api_calls']:
    print(f"  {call['timestamp']}: {call['input_tokens']} in, {call['output_tokens']} out")
```

### Aggregate Results

```python
# Load all annotated datasets from a run
all_datasets = []

for tissue_dir in run_dir.iterdir():
    if tissue_dir.is_dir() and tissue_dir.name not in ['batch_summary.txt', 'batch_summary.csv']:
        for dataset_dir in tissue_dir.iterdir():
            if dataset_dir.is_dir():
                h5ad_file = dataset_dir / f"{dataset_dir.name}_annotated.h5ad"
                if h5ad_file.exists():
                    adata = sc.read_h5ad(h5ad_file)
                    adata.obs['dataset'] = dataset_dir.name
                    adata.obs['tissue'] = tissue_dir.name
                    all_datasets.append(adata)

# Concatenate all datasets
combined = sc.concat(all_datasets, join='outer')
print(f"Combined: {combined.n_obs:,} cells from {len(all_datasets)} datasets")
```

### Compare Runs

```bash
# List all runs (organized by tissue)
ls -lh batch_results/

# List only brain runs
ls -d batch_results/run_brain_*/

# Compare token usage across runs
for run in batch_results/run_*/; do
    echo "=== $run ==="
    grep "Total cost:" "$run/batch_summary.txt"
done

# Compare number of files processed
for run in batch_results/run_*/; do
    files=$(grep -c "success" "$run/batch_summary.csv" || echo "0")
    echo "$run: $files files"
done

# Find all brain processing runs
ls -d batch_results/run_brain_* 2>/dev/null || echo "No brain runs found"
```

## Tips

1. **Archive old runs** - Move or compress old `run_*` directories to save space
   ```bash
   # Archive all brain runs
   tar -czf archive_brain_runs.tar.gz batch_results/run_brain_*/
   
   # Archive by date
   tar -czf archive_20260216.tar.gz batch_results/run_*_20260216_*/
   ```

2. **Quick stats** - Use `batch_summary.txt` for quick overview
   ```bash
   cat batch_results/run_brain_20260216_151425/batch_summary.txt
   ```

3. **Find specific datasets** - Use `find` to locate files
   ```bash
   find batch_results/run_brain_20260216_151425 -name "*Bassez*_annotated.h5ad"
   ```

4. **Filter runs by tissue** - Use wildcards to find specific tissue runs
   ```bash
   # Show only brain runs
   ls -lh batch_results/run_brain_*/
   
   # Show runs with multiple tissues
   ls -lh batch_results/run_*_*_2026*/
   ```

4. **Check for failures** - Filter the CSV for failed files
   ```python
   failed = summary[summary['status'] == 'failed']
   print(failed[['dataset', 'tissue', 'error']])
   ```

## Best Practices

- **Don't modify run directories** - They represent a complete, immutable record of a processing run
- **Use descriptive output-dir names** for different experiments:
  ```bash
  python batch_process_all_tissues.py --output-dir batch_results_pilot_study
  ```
- **Keep token_tracking.json** for cost analysis and budgeting
- **Document your runs** - Add notes about parameters used:
  ```bash
  echo "Testing resolution=0.8 with 50 PCs" > batch_results/run_20260216_151425/NOTES.txt
  ```
