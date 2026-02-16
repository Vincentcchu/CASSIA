# CASSIA Batch Processing

This directory contains a script for batch processing multiple h5ad files with CASSIA cell type annotation.

## Overview

The `batch_process_all_tissues.py` script provides a flexible CLI for processing:
- All tissues in your data directory
- A single specific tissue
- Multiple specific tissues
- Specific files matching a pattern

## Requirements

- CASSIA installed (with token tracking support)
- scanpy, pandas, numpy installed
- API keys configured in `../api_keys.json`

## Quick Start

### Preview Files Without Processing

```bash
# Preview all files
python batch_process_all_tissues.py --preview

# Preview a specific tissue
python batch_process_all_tissues.py --tissue brain --preview

# Preview with file filter
python batch_process_all_tissues.py --tissue breast --files "Data_Bassez*.h5ad" --preview
```

### Process Files

```bash
# Process all tissues
python batch_process_all_tissues.py

# Process a single tissue
python batch_process_all_tissues.py --tissue brain

# Process multiple specific tissues
python batch_process_all_tissues.py --tissues brain breast colorectal

# Process specific files in a tissue
python batch_process_all_tissues.py --tissue breast --files "Data_Bassez*.h5ad"
```

## Command-Line Options

### Required Options (choose one)

```bash
--tissue TISSUE               # Process a single tissue (e.g., brain, breast)
--tissues TISSUE1 TISSUE2 ... # Process multiple tissues
# (no flag)                   # Process all tissues in data directory
```

### File Selection

```bash
--files "PATTERN"             # File pattern to match (default: *.h5ad)
                              # Examples: "Data_Bassez*.h5ad", "Data_*.h5ad"
```

### Data Location

```bash
--data-root PATH              # Root data directory
                              # Default: /cs/student/projects2/aisd/2024/shekchu/projects/data
```

### Model Configuration

```bash
--provider {openai,anthropic,openrouter}  # LLM provider (default: openai)
--model MODEL                             # Model name (default: gpt-4o-2024-08-06)
--species SPECIES                         # Species (default: Human)
```

### Processing Parameters

```bash
--resolution FLOAT            # Clustering resolution (default: 0.5)
--n-pcs INT                   # Number of principal components (default: 40)
--n-neighbors INT             # Number of neighbors (default: 10)
--max-workers INT             # Maximum parallel workers (default: 4)
--score-threshold FLOAT       # Score threshold for boost annotation (default: 75)
```

### Output Configuration

```bash
--output-dir PATH             # Output directory (default: ./batch_results)
```

### Cost Tracking

```bash
--input-cost FLOAT            # Input cost per 1M tokens in USD (default: 1.25)
--output-cost FLOAT           # Output cost per 1M tokens in USD (default: 10.0)
```

### Preview Mode

```bash
--preview                     # Preview files without processing
                              # Shows: files found, cell/gene counts, time/cost estimates
```

## Examples

### Example 1: Preview All Brain Files

```bash
python batch_process_all_tissues.py --tissue brain --preview
```

Output:
```
======================================================================
PREVIEW MODE - Files to be processed:
======================================================================

BRAIN (8 files):
   1. Data_Choudhury2022_Brain.h5ad
      5,432 cells × 20,123 genes
   2. Data_Filbin2018_Brain.h5ad
      2,587 cells × 23,686 genes
   ...

======================================================================
ESTIMATES:
======================================================================
Time: ~64 minutes (1.1 hours)
Tokens: ~1.2M tokens
Cost: ~$5.23
```

### Example 2: Process Specific Files

```bash
python batch_process_all_tissues.py --tissue breast --files "Data_Bassez*.h5ad"
```

### Example 3: Process Multiple Tissues

```bash
python batch_process_all_tissues.py --tissues brain breast colorectal
```

### Example 4: Custom Model Configuration

```bash
python batch_process_all_tissues.py \
    --tissue lung \
    --provider openrouter \
    --model "google/gemini-2.5-flash" \
    --resolution 0.8 \
    --n-pcs 50
```

## Output Structure

Results are organized in tissue-labeled timestamped run directories to keep each batch processing session separate:

```
batch_results/
└── run_brain_20260216_151425/                   # Tissue + timestamp run directory
    ├── batch_summary.txt                        # Human-readable summary report
    ├── batch_summary.csv                        # Detailed results table (for analysis)
    ├── token_tracking.json                      # Token usage and API call details
    ├── brain/                                   # Tissue-specific folder
    │   ├── Data_Choudhury2022_Brain/            # Dataset-specific folder
    │   │   ├── Data_Choudhury2022_Brain_markers.csv
    │   │   ├── Data_Choudhury2022_Brain_annotated.h5ad
    │   │   ├── Data_Choudhury2022_Brain_annotation.txt
    │   │   └── Data_Choudhury2022_Brain_summary.txt
    │   ├── Data_Filbin2018_Brain/
    │   │   └── ...
    │   └── ...
    └── breast/
        ├── Data_Bassez2021_Breast/
        │   └── ...
        └── ...
```

### Output Files Explained

**Run-level files** (in `run_YYYYMMDD_HHMMSS/`):
- `batch_summary.txt` - Overall summary with timing, token usage, and costs
- `batch_summary.csv` - Structured data for all processed files (for pandas/Excel)
- `token_tracking.json` - Detailed token tracking (input/output per API call)

**Dataset-level files** (in `run_YYYYMMDD_HHMMSS/{tissue}/{dataset}/`):
- `{dataset}_markers.csv` - Marker genes for each cluster (DEG analysis results)
- `{dataset}_annotated.h5ad` - Annotated AnnData object with cell type labels
- `{dataset}_annotation.txt` - CASSIA's text-format annotation output
- `{dataset}_summary.txt` - Individual file processing summary (if generated)

### Naming Convention

- **Dataset names** are extracted from h5ad filenames
  - Example: `Data_Bassez2021_Breast.h5ad` → dataset name is `Data_Bassez2021_Breast`
- **All output files** for a dataset use this consistent prefix
- **Run directories** use timestamp format `run_YYYYMMDD_HHMMSS` to:
  - Avoid overwriting previous runs
  - Make it easy to track when processing occurred
  - Allow multiple runs with different parameters

### Locating Your Results

After running the script, all outputs are in a single timestamped directory. The script will print:
```
Output directory: batch_results/run_20260216_151425
```

You can then:
```bash
# View the summary
cat batch_results/run_20260216_151425/batch_summary.txt

# Load results in Python
import pandas as pd
results = pd.read_csv('batch_results/run_20260216_151425/batch_summary.csv')

# Access annotated data
import scanpy as sc
adata = sc.read_h5ad('batch_results/run_20260216_151425/brain/Data_Filbin2018_Brain/Data_Filbin2018_Brain_annotated.h5ad')
```

## Summary Report

Each run generates a summary showing:
- Total files processed
- Total time elapsed
- Token usage (input/output/total)
- Total cost
- Per-file breakdowns

## Expected Directory Structure

The script expects your data to be organized as:

```
data/
├── brain/
│   └── h5ad_unlabelled/
│       ├── Data_Choudhury2022_Brain.h5ad
│       ├── Data_Filbin2018_Brain.h5ad
│       └── ...
├── breast/
│   └── h5ad_unlabelled/
│       ├── Data_Bassez2021_Breast.h5ad
│       └── ...
└── colorectal/
    └── h5ad_unlabelled/
        └── ...
```

## Tips

1. **Always preview first** to verify files and estimate costs:
   ```bash
   python batch_process_all_tissues.py --preview
   ```

2. **Test with one file** before batch processing:
   ```bash
   python batch_process_all_tissues.py --tissue brain --files "Data_Filbin*.h5ad"
   ```

3. **Monitor memory usage** for large datasets - adjust `--max-workers` if needed

4. **Check the summary file** after completion for detailed statistics

## Troubleshooting

### No files found

- Verify your `--data-root` path is correct
- Check that `h5ad_unlabelled/` subdirectories exist
- Verify file pattern with `--preview`

### API errors

- Ensure `../api_keys.json` exists and has correct keys
- Check provider and model names are correct

### Memory issues

- Reduce `--max-workers` (default: 4)
- Process one tissue at a time

## Notes

- The script automatically tracks token usage and calculates costs
- Processing time depends on dataset size and model speed
- All timing and token metrics are saved in summary files
- Failed files are logged but don't stop batch processing
