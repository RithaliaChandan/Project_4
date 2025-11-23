Project — TCGA-LUAD NKX2-1 Expression Analysis
Author: Chandan Rithalia
Date:23/11/2025

Objective:
To analyze NKX2-1 (TTF-1) gene expression differences between Primary Tumor and Solid Tissue Normal lung samples in TCGA-LUAD, using TPM values derived from RNA-seq augmented STAR gene count files. The goal is to visualize expression variation by generating a tumor vs. normal boxplot.

Steps:
1. Set working directory and extract tcga_data.tar.gz.
2. Load gdc_sample_sheet.tsv and filter tumor/normal samples.
3. Keep only RNA-seq augmented STAR gene count files.
4. Extract NKX2-1 TPM values from each sample folder.
5. Compute log2(TPM + 1).
6. Remove missing TPM values.
7. Plot tumor vs. normal NKX2-1 expression.

Scripts Used:
- get_tpm_NKX2 function (extracts TPM values)
- Main analysis R script (data filtering, TPM extraction, plot generation)

Input Files:
- tcga_data.tar.gz — unpacked RNA-seq data folders
- gdc_sample_sheet.tsv — sample metadata

Output Files:
- NKX2-1_boxplot.png (or R plot window output)
- ss_clean dataframe with TPM values
- log2(TPM + 1) column

Conclusion:
The resulting boxplot confirms that NKX2-1 expression is significantly higher in tumor samples compared to normal lung tissues, supporting its known role as a diagnostic and functional biomarker in lung adenocarcinoma.
