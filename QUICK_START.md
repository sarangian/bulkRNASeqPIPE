# Quick Start Guide: Enhanced DESeq2 Script

## Installation

```r
# Install the new required package
install.packages("openxlsx")
```

That's it! All other dependencies remain the same.

---

## Running the Script

**No changes to command line usage!**

```bash
Rscript DESeq2_enhanced.R \
  --projectName MyAnalysis \
  --design design.txt \
  --featureQuant featureCounts.txt \
  --varInt group \
  --condControl control \
  --condTreated treated
```

---

## What You Get (NEW)

### 1. Enhanced Results File
**Location:** `tables/*_complete_with_counts.tsv`

**Contains:**
- Raw counts for all samples
- Normalized counts for all samples  
- Complete DESeq2 statistics

### 2. Six Filtered Result Files  
**Location:** `filtered_results/`

**Files:**
- `*_unadjusted_p0.05_lfc2.tsv` (most stringent, exploratory)
- `*_unadjusted_p0.05_lfc1.5.tsv`
- `*_unadjusted_p0.05_lfc1.tsv`
- `*_adjusted_p0.05_lfc2.tsv` (**recommended for publication**)
- `*_adjusted_p0.05_lfc1.5.tsv` (**recommended for publication**)
- `*_adjusted_p0.05_lfc1.tsv` (**recommended for publication, most sensitive**)

### 3. Excel Workbook
**Location:** `filtered_results/*_all_cutoffs.xlsx`

**Features:**
- All 6 cutoffs as separate tabs
- Formatted headers
- Ready to share

### 4. Summary Table
**Location:** `tables/*_cutoff_summary.txt`

**Shows:** DEG counts for all 6 cutoffs

---

## Which Cutoff Should I Use?

### For Publication 📄
**Use:** `adjusted_p0.05_lfc1.tsv` or `adjusted_p0.05_lfc1.5.tsv`
- Adjusted p-values control false discovery rate
- Standard |log2FC| ≥ 1 (2-fold change)
- Increase to 1.5 if you have too many DEGs

### For qRT-PCR Validation 🔬
**Use:** `adjusted_p0.05_lfc2.tsv`
- High fold-change genes
- More reliable for validation
- Fewer genes to validate

### For Exploratory Analysis 🔍
**Use:** `unadjusted_p0.05_lfc1.5.tsv` or `unadjusted_p0.05_lfc2.tsv`
- Cast a wider net
- Generate hypotheses
- More false positives expected

### For Pathway Analysis 🗺️
**Use:** `adjusted_p0.05_lfc1.tsv`
- More genes = better pathway coverage
- Still controlled FDR
- Standard threshold

---

## Quick Checks After Running

### 1. Check the Cutoff Summary
```bash
cat MyAnalysis/tables/*_cutoff_summary.txt
```

**Look for:**
- How many DEGs in each cutoff?
- Big difference between adjusted and unadjusted?
- Are the numbers reasonable for your experiment?

### 2. Open the Excel File
```bash
# Open with Excel, LibreOffice, or:
open MyAnalysis/filtered_results/*_all_cutoffs.xlsx
```

**Compare:**
- Switch between tabs to see different cutoffs
- Notice how gene lists change with stringency
- Pick the most appropriate for your needs

### 3. Check Sample Counts in Results
```bash
head MyAnalysis/tables/*_complete_with_counts.tsv
```

**Verify:**
- Raw counts look reasonable?
- Normalized counts similar across replicates?
- Any obvious outliers?

---

## Common Workflows

### Workflow 1: Standard Analysis

```bash
# 1. Run analysis
Rscript DESeq2_enhanced.R \
  --design design.txt \
  --featureQuant counts.txt \
  --varInt group \
  --condControl control \
  --condTreated treated

# 2. Check summary
cat MyProject/tables/*_cutoff_summary.txt

# 3. Choose appropriate cutoff (e.g., adjusted p, |LFC| ≥ 1.5)
# 4. Use that TSV file for downstream analysis
```

### Workflow 2: Share with Collaborators

```bash
# 1. Run analysis (same as above)

# 2. Send them the Excel file
# File: MyProject/filtered_results/*_all_cutoffs.xlsx

# 3. They can:
#    - Switch between tabs for different stringencies
#    - Filter/sort in Excel
#    - Export genes they're interested in
```

### Workflow 3: Pathway Enrichment

```bash
# 1. Run analysis

# 2. Extract gene list from appropriate cutoff
cut -f1 MyProject/filtered_results/*_adjusted_p0.05_lfc1.tsv | tail -n +2 > gene_list.txt

# 3. Upload gene_list.txt to:
#    - DAVID (https://david.ncifcrf.gov/)
#    - Enrichr (https://maayanlab.cloud/Enrichr/)
#    - g:Profiler (https://biit.cs.ut.ee/gprofiler/)
```

---

## Understanding the Results

### TSV File Structure

```
Gene | raw_Sample1 | raw_Sample2 | ... | norm_Sample1 | norm_Sample2 | ... | 
baseMean | log2FoldChange | pvalue | padj | FoldChange | Significance | Regulation
```

**Key columns:**
- `Gene`: Gene identifier
- `raw_*`: Raw counts per sample
- `norm_*`: Normalized counts per sample
- `log2FoldChange`: Effect size (+ = upregulated, - = downregulated)
- `pvalue`: Nominal p-value
- `padj`: Adjusted p-value (FDR)
- `FoldChange`: Linear fold change (2^log2FC)
- `Regulation`: Up-regulated / Down-regulated

### Excel Worksheets

Each worksheet contains:
1. Same structure as TSV files
2. Formatted header (blue, bold, white text)
3. Frozen first row
4. All DEGs passing that specific cutoff

---

## Tips

### Tip 1: Start with Summary
Always check `*_cutoff_summary.txt` first to get an overview of DEG counts.

### Tip 2: Use Adjusted P-values for Final Results
Unadjusted p-values are for exploration only. Publications should use adjusted p-values (FDR).

### Tip 3: Consider Biological Context
- Transcription factors often have smaller fold changes
- Structural genes may show larger fold changes
- Adjust |log2FC| threshold accordingly

### Tip 4: Excel for Quick Review
Open the Excel file to quickly browse through different cutoffs without re-running analysis.

### Tip 5: Keep Raw Data
The enhanced results file with sample counts is great for:
- Troubleshooting outliers
- Making custom plots
- Checking specific genes across samples

---

## Example Output Structure

```
MyAnalysis/
├── tables/
│   ├── treated_vs_control_complete_with_counts.tsv    ← All data
│   ├── treated_vs_control_cutoff_summary.txt          ← Quick overview
│   └── ... (other files)
│
└── filtered_results/
    ├── treated_vs_control_unadjusted_p0.05_lfc2.tsv   ← Individual cutoffs
    ├── treated_vs_control_unadjusted_p0.05_lfc1.5.tsv
    ├── treated_vs_control_unadjusted_p0.05_lfc1.tsv
    ├── treated_vs_control_adjusted_p0.05_lfc2.tsv     ← Recommended
    ├── treated_vs_control_adjusted_p0.05_lfc1.5.tsv   ← Recommended
    ├── treated_vs_control_adjusted_p0.05_lfc1.tsv     ← Recommended
    └── treated_vs_control_all_cutoffs.xlsx            ← Excel workbook
```

---

## Troubleshooting

### "Package 'openxlsx' not found"
```r
install.packages("openxlsx")
```

### Too Many DEGs
- Use more stringent cutoff (higher |log2FC|)
- Use adjusted p-values instead of unadjusted
- Check if samples are correctly assigned to groups

### Too Few DEGs
- Check if biological replicates are sufficient (≥3 per group)
- Check if sequencing depth is adequate
- Consider using less stringent cutoff for exploration
- Use unadjusted p-values with high |log2FC| threshold

### Excel File Won't Open
- Try LibreOffice Calc instead of Microsoft Excel
- Check file wasn't corrupted during transfer
- Ensure filename/path isn't too long

---

## Command Line Options (Quick Reference)

All original options work the same:

```bash
--projectName      # Project name (default: current directory)
--design           # Design/target file
--featureQuant     # FeatureCounts output
--varInt           # Variable of interest (default: group)
--condControl      # Control condition
--condTreated      # Treatment condition
--batch            # Batch variable (optional)
--filterMethod     # HTSFilter/genefilter/none (default: HTSFilter)
--alpha            # FDR threshold (default: 0.05)
--pAdjustMethod    # BH/BY (default: BH)
--fitType          # parametric/local/mean (default: parametric)
--ntop             # Top genes for PCA (default: 500)
```

---

## What Didn't Change

✅ All command-line options (same syntax)
✅ All QC plots  
✅ All statistical methods
✅ Gene filtering options
✅ 3D interactive PCA
✅ Original result files still generated
✅ Output directory structure (just added `filtered_results/`)

**The enhanced script is 100% backward compatible!**

---

## Next Steps

1. **Run the script** with your data
2. **Check cutoff_summary.txt** to see DEG counts
3. **Open the Excel file** to browse results
4. **Choose appropriate cutoff** for your analysis
5. **Use TSV files** for downstream analysis

---

## Support

For questions or issues:
- Check ENHANCEMENT_DOCUMENTATION.md for details
- Review TROUBLESHOOTING.md for common issues
- Contact: BaseSolve Informatics Pvt. Ltd, Ahmedabad, Gujarat

---

## Summary

**What you need to remember:**

1. ✅ Install `openxlsx` package
2. ✅ Run script exactly as before (no changes)
3. ✅ You get 6 filtered result files automatically
4. ✅ Use adjusted p-value cutoffs for publication
5. ✅ Excel file has all cutoffs for easy comparison

**That's it! Enjoy the enhanced analysis!** 🎉
