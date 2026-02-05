# DESeq2 Enhanced Script Documentation

## Overview

This enhanced version of the DESeq2 script maintains **ALL original functionalities** while adding powerful new features for comprehensive differential expression analysis.

---

## What's New? 🆕

### 1. **Enhanced Results with Sample Counts**

**Original:** Results contained only DESeq2 statistics (log2FC, p-values, etc.)

**Enhanced:** Results now include:
- ✅ Raw counts for **every sample**
- ✅ Normalized counts for **every sample**
- ✅ All DESeq2 statistics
- ✅ Fold change (2^log2FC)
- ✅ Significance classification
- ✅ Regulation direction

**Output File:** `{comparison}_complete_with_counts.tsv`

**Example columns:**
```
Gene | raw_Sample1 | raw_Sample2 | ... | norm_Sample1 | norm_Sample2 | ... | 
baseMean | log2FoldChange | pvalue | padj | FoldChange | Significance | Regulation
```

---

### 2. **Six Automatic Filtering Cutoffs**

The script now automatically applies **6 different combinations** of p-value and fold-change thresholds:

| # | P-value Type | Threshold | Log2FC Threshold | Use Case |
|---|--------------|-----------|------------------|----------|
| 1 | Unadjusted | ≤ 0.05 | ≥ 2.0 | High stringency, exploratory |
| 2 | Unadjusted | ≤ 0.05 | ≥ 1.5 | Moderate stringency, exploratory |
| 3 | Unadjusted | ≤ 0.05 | ≥ 1.0 | Lower stringency, exploratory |
| 4 | **Adjusted (FDR)** | ≤ 0.05 | ≥ 2.0 | **Publication (high effect)** |
| 5 | **Adjusted (FDR)** | ≤ 0.05 | ≥ 1.5 | **Publication (moderate effect)** |
| 6 | **Adjusted (FDR)** | ≤ 0.05 | ≥ 1.0 | **Publication (standard)** |

**Recommended for publications:** Cutoffs 4-6 (adjusted p-values control false discovery rate)

---

### 3. **TSV Output for Each Cutoff**

**What:** Individual tab-separated value files for each cutoff

**Location:** `{projectName}/filtered_results/`

**Files created:**
```
treated_vs_control_unadjusted_p0.05_lfc2.tsv
treated_vs_control_unadjusted_p0.05_lfc1.5.tsv
treated_vs_control_unadjusted_p0.05_lfc1.tsv
treated_vs_control_adjusted_p0.05_lfc2.tsv
treated_vs_control_adjusted_p0.05_lfc1.5.tsv
treated_vs_control_adjusted_p0.05_lfc1.tsv
```

**Use:** Easy import into:
- R for further analysis
- Python (pandas)
- GraphPad Prism
- Pathway analysis tools (DAVID, Enrichr, etc.)

---

### 4. **Excel Workbook with All Cutoffs**

**What:** Single Excel file with 6 worksheets (one per cutoff)

**Location:** `{projectName}/filtered_results/{comparison}_all_cutoffs.xlsx`

**Features:**
- ✅ Professional formatted headers (blue background, white text, bold)
- ✅ Frozen first row for easy scrolling
- ✅ Auto-adjusted column widths
- ✅ All 6 cutoffs in one file for easy comparison
- ✅ Ready to share with collaborators

**Worksheets:**
1. `unadjusted_p0.05_lfc2`
2. `unadjusted_p0.05_lfc1.5`
3. `unadjusted_p0.05_lfc1`
4. `adjusted_p0.05_lfc2`
5. `adjusted_p0.05_lfc1.5`
6. `adjusted_p0.05_lfc1`

---

### 5. **Comprehensive Cutoff Summary**

**What:** Summary table showing DEG counts for all cutoffs

**Location:** `{projectName}/tables/{comparison}_cutoff_summary.txt`

**Example:**
```
Cutoff_Name              Description                                   Total_DEGs  Upregulated  Downregulated
unadjusted_p0.05_lfc2   Unadjusted p-value ≤ 0.05, |log2FC| ≥ 2      150         85           65
unadjusted_p0.05_lfc1.5 Unadjusted p-value ≤ 0.05, |log2FC| ≥ 1.5    234         128          106
unadjusted_p0.05_lfc1   Unadjusted p-value ≤ 0.05, |log2FC| ≥ 1      456         245          211
adjusted_p0.05_lfc2     Adjusted p-value ≤ 0.05, |log2FC| ≥ 2        89          48           41
adjusted_p0.05_lfc1.5   Adjusted p-value ≤ 0.05, |log2FC| ≥ 1.5      145         78           67
adjusted_p0.05_lfc1     Adjusted p-value ≤ 0.05, |log2FC| ≥ 1        298         161          137
```

**Use:** Quick overview to choose appropriate cutoff for your study

---

## Original Features Retained ✅

All original functionalities are **preserved and unchanged**:

### Gene Filtering
- ✅ HTSFilter (Jaccard similarity-based)
- ✅ genefilter (iterative variance-based)
- ✅ No filtering option (DESeq2 independent filtering only)

### Statistical Analysis
- ✅ Full DESeq2 pipeline with proper normalization
- ✅ Dispersion estimation (parametric, local, or mean)
- ✅ Multiple testing correction (BH or BY)
- ✅ Batch effect correction support

### Quality Control
- ✅ Sample distance heatmap
- ✅ 2D PCA plots (PNG and SVG)
- ✅ **3D interactive PCA** (plotly)
- ✅ Dispersion plots
- ✅ Cook's distance QC

### Visualizations
- ✅ Enhanced volcano plots (EnhancedVolcano)
- ✅ MA plots (DESeq2 native)
- ✅ VIDGER plots (MA, volcano, box, scatter)
- ✅ Gene heatmaps (top DEGs)
- ✅ Interactive Glimma plots (volcano, MDS, MA)

### Results Export
- ✅ Complete DESeq2 results
- ✅ Significant genes only
- ✅ Filtered by log2FC thresholds
- ✅ PCA coordinates
- ✅ Session information

---

## Installation

### New Dependency

The enhanced script requires one additional package:

```r
install.packages("openxlsx")
```

### Complete Installation

```r
# Install all required packages
install.packages(c(
  "optparse", "ggplot2", "dplyr", "readr",
  "RColorBrewer", "pheatmap", "reshape",
  "ggrepel", "plotly", "htmlwidgets",
  "matrixStats", "openxlsx"  # NEW
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "HTSFilter", "genefilter",
  "vidger", "EnhancedVolcano",
  "GenomicFeatures", "regionReport",
  "DEFormats", "Glimma", "DelayedArray"
))
```

---

## Usage

### Basic Usage (Same as Before)

```bash
Rscript DESeq2_enhanced.R \
  --projectName MyRNAseq \
  --design design.txt \
  --featureQuant featureCounts.txt \
  --varInt group \
  --condControl control \
  --condTreated treated \
  --filterMethod HTSFilter \
  --alpha 0.05
```

**No changes needed to command line arguments!**

---

## Output Files Comparison

### Original Script Output
```
MyProject/
├── tables/
│   ├── treated_vs_control_complete_DESeq2_results.tsv
│   ├── treated_vs_control_significant_genes_padj0.05.tsv
│   ├── gene_filtering_summary.txt
│   ├── PCA_coordinates.tsv
│   └── analysis_summary.txt
├── figures/
│   ├── enhanced_volcano_plot.png
│   ├── MA_plot_DESeq2.png
│   └── (other plots)
├── QC/
│   └── (QC plots)
└── html/
    └── (interactive plots)
```

### Enhanced Script Output (Additional Files)
```
MyProject/
├── tables/
│   ├── ... (all original files) ...
│   ├── treated_vs_control_complete_with_counts.tsv      ← NEW
│   └── treated_vs_control_cutoff_summary.txt            ← NEW
└── filtered_results/                                     ← NEW DIRECTORY
    ├── treated_vs_control_unadjusted_p0.05_lfc2.tsv
    ├── treated_vs_control_unadjusted_p0.05_lfc1.5.tsv
    ├── treated_vs_control_unadjusted_p0.05_lfc1.tsv
    ├── treated_vs_control_adjusted_p0.05_lfc2.tsv
    ├── treated_vs_control_adjusted_p0.05_lfc1.5.tsv
    ├── treated_vs_control_adjusted_p0.05_lfc1.tsv
    └── treated_vs_control_all_cutoffs.xlsx              ← Excel workbook
```

---

## Workflow Recommendations

### 1. Standard Analysis Workflow

```bash
# Run the enhanced script
Rscript DESeq2_enhanced.R \
  --design design.txt \
  --featureQuant counts.txt \
  --varInt group \
  --condControl control \
  --condTreated treated
```

### 2. Review Results

1. **Check cutoff summary** (`tables/*_cutoff_summary.txt`)
   - Compare DEG counts across all cutoffs
   - Identify appropriate stringency for your study

2. **Open Excel file** (`filtered_results/*_all_cutoffs.xlsx`)
   - Quick visual review of DEGs at different cutoffs
   - Easy to share with collaborators

3. **Choose appropriate cutoff**
   - **For publication:** Use adjusted p-value cutoffs (4-6)
   - **For exploration:** Use unadjusted p-value cutoffs (1-3)
   - **For biological significance:** Consider log2FC threshold

### 3. Downstream Analysis

**Option A:** Use TSV files directly
```r
# Read specific cutoff results
degs <- read.table("filtered_results/treated_vs_control_adjusted_p0.05_lfc1.tsv", 
                   header = TRUE, sep = "\t")

# Extract gene list for pathway analysis
gene_list <- degs$Gene
```

**Option B:** Use complete results with counts
```r
# Read enhanced results
results <- read.table("tables/treated_vs_control_complete_with_counts.tsv",
                      header = TRUE, sep = "\t")

# Access raw counts
raw_counts <- results[, grep("^raw_", colnames(results))]

# Access normalized counts
norm_counts <- results[, grep("^norm_", colnames(results))]
```

---

## Understanding the Cutoffs

### Unadjusted vs. Adjusted P-values

**Unadjusted P-values (cutoffs 1-3):**
- Nominal p-values without multiple testing correction
- Higher sensitivity but more false positives
- Use for: Exploratory analysis, hypothesis generation
- ⚠️ Risk: ~5% false positive rate per test

**Adjusted P-values / FDR (cutoffs 4-6):**
- Corrected for multiple testing (Benjamini-Hochberg by default)
- Controls false discovery rate across all tests
- Use for: Publication-quality results, validation experiments
- ✅ Recommended for final gene lists

### Log2 Fold Change Thresholds

| Threshold | Linear FC | Interpretation |
|-----------|-----------|----------------|
| ±1.0 | 2-fold | Standard threshold, biologically relevant |
| ±1.5 | 2.83-fold | Moderate effect size |
| ±2.0 | 4-fold | Large effect size, highly stringent |

**Recommendation:**
- Start with |log2FC| ≥ 1 for broad overview
- Increase to 1.5 or 2 for more robust gene lists
- Consider biological context (e.g., transcription factors may have smaller changes)

---

## Practical Examples

### Example 1: Standard RNA-seq Analysis

**Goal:** Identify DEGs between treatment and control

```bash
Rscript DESeq2_enhanced.R \
  --design design.txt \
  --featureQuant counts.txt \
  --varInt condition \
  --condControl untreated \
  --condTreated drug_treatment \
  --filterMethod HTSFilter \
  --alpha 0.05
```

**Recommended cutoff:** `adjusted_p0.05_lfc1.tsv` (standard for publication)

---

### Example 2: Time-Series Experiment

**Goal:** Compare multiple timepoints

```bash
# Note: For multiple conditions, run pairwise comparisons
# Example: Compare each timepoint to baseline

# Timepoint 1 vs Baseline
Rscript DESeq2_enhanced.R \
  --design design_t1.txt \
  --condControl baseline \
  --condTreated timepoint1 \
  ...

# Timepoint 2 vs Baseline
Rscript DESeq2_enhanced.R \
  --design design_t2.txt \
  --condControl baseline \
  --condTreated timepoint2 \
  ...
```

**Recommended cutoff:** `adjusted_p0.05_lfc1.5.tsv` (moderate stringency)

---

### Example 3: Exploratory Analysis with Many Samples

**Goal:** Find candidate genes for validation

```bash
Rscript DESeq2_enhanced.R \
  --design design.txt \
  --featureQuant counts.txt \
  --filterMethod genefilter \
  --genefilterRetention 0.50 \
  --alpha 0.1  # More relaxed
  ...
```

**Recommended cutoff:** `unadjusted_p0.05_lfc2.tsv` (high fold-change, exploratory)

---

## Tips and Best Practices

### 1. Choosing the Right Cutoff

✅ **For publications:**
- Use adjusted p-values (FDR)
- Start with |log2FC| ≥ 1
- Check how many DEGs you get
- Increase to 1.5 or 2 if list is too large

✅ **For qRT-PCR validation:**
- Use adjusted p-values
- Consider |log2FC| ≥ 1.5 or 2
- Select genes with high effect sizes

✅ **For pathway analysis:**
- Use adjusted p-values
- Can use |log2FC| ≥ 1
- More genes = better pathway coverage

### 2. Working with Excel Files

The Excel workbook allows easy comparison:
1. Open `*_all_cutoffs.xlsx`
2. Switch between tabs to see different cutoffs
3. Use Excel's built-in filtering and sorting
4. Copy genes to separate sheets for validation lists

### 3. Extracting Specific Gene Lists

```r
# Load a specific cutoff
degs <- read.table("filtered_results/treated_vs_control_adjusted_p0.05_lfc1.5.tsv",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Get upregulated genes only
up_genes <- degs[degs$Regulation == "Up-regulated", "Gene"]

# Get top 20 by fold change
top20 <- head(degs[order(-abs(degs$log2FoldChange)), ], 20)

# Export for pathway analysis
write.table(degs$Gene, "gene_list_for_enrichment.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
```

### 4. Comparing Cutoffs

Use the summary file to decide:

```r
summary <- read.table("tables/treated_vs_control_cutoff_summary.txt",
                      header = TRUE, sep = "\t")

# Visualize DEG counts
barplot(summary$Total_DEGs, names.arg = 1:6,
        main = "DEGs Across Cutoffs",
        xlab = "Cutoff", ylab = "Number of DEGs")
```

---

## Troubleshooting

### Issue: Too Many DEGs in Unadjusted Results

**Solution:** Use adjusted p-values (cutoffs 4-6) or increase log2FC threshold

### Issue: Too Few DEGs in Adjusted Results

**Check:**
1. Sample size adequate? (Need ≥3 replicates per group)
2. Biological variation high?
3. Sequencing depth sufficient?

**Consider:**
- Using unadjusted p-values for exploratory analysis
- Relaxing alpha to 0.10
- Focusing on high fold-change genes

### Issue: Excel File Won't Open

**Possible causes:**
1. File path too long → Move to shorter path
2. File corrupted → Re-run analysis
3. Excel version incompatible → Try LibreOffice Calc

---

## Performance Notes

**Memory usage:**
- Enhanced results use more memory (includes all sample counts)
- For very large datasets (>50,000 genes, >100 samples), consider:
  - Using filtered gene sets
  - Running on machine with ≥16GB RAM

**Runtime:**
- Additional processing adds ~10-30 seconds
- Excel creation adds ~5-10 seconds
- Total impact: <1 minute for typical datasets

---

## Comparison with Original Script

| Feature | Original | Enhanced |
|---------|----------|----------|
| Core DESeq2 analysis | ✅ | ✅ |
| Gene filtering options | ✅ | ✅ |
| QC plots | ✅ | ✅ |
| 3D PCA | ✅ | ✅ |
| Results with sample counts | ❌ | ✅ |
| Multiple cutoff filtering | ❌ | ✅ |
| TSV export per cutoff | ❌ | ✅ |
| Excel workbook | ❌ | ✅ |
| Cutoff summary table | ❌ | ✅ |

---

## Version Information

**Enhanced Version:** 2.0
**Based on:** DESeq2 pipeline v1.0
**Date:** February 2026
**Compatibility:** Fully backward compatible

---

## Support

**BaseSolve Informatics Pvt. Ltd**
Ahmedabad, Gujarat, India

For questions about:
- Script usage and parameters
- Interpreting results
- Troubleshooting errors
- Custom analysis requests

---

## Summary

The enhanced script provides:

1. ✅ **All original features** preserved
2. ✅ **Sample-level counts** in results
3. ✅ **Six automatic cutoffs** for flexible analysis
4. ✅ **Publication-ready Excel files**
5. ✅ **Easy-to-use TSV files** for downstream analysis
6. ✅ **No changes** to command-line usage

**Recommendation:** Use this enhanced version for all new analyses. It provides more comprehensive results with minimal additional runtime.
