#Set working directory
setwd("C:/Users/Chandan/Downloads/TCGA_LUAD")

untar("tcga_data.tar.gz")
ss <- read.delim("gdc_sample_sheet.tsv", stringsAsFactors = FALSE)

ss_sub <- subset(ss,Sample.Type %in% c("Primary Tumor", "Solid Tissue Normal"))
ss_sub$Dir <- ss_sub$File.ID

# Keep only RNA-seq augmented STAR gene count files
rna_rows <- grepl("rna_seq\\.augmented_star_gene_counts\\.tsv$", ss_sub$File.Name)
ss_rna <- ss_sub[rna_rows, ]

#extract NKX2-1 TPM from a single RNA-seq folder
get_tpm_NKX2 <- function(folder_name) {
  file_path <- list.files(
    folder_name,
    pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$",
    full.names = TRUE )
  if (length(file_path) == 0) return(NA_real_)
  
# Read the raw lines
  lines <- readLines(file_path[1])
  header_idx <- grep("^gene_id", lines)[1]
  if (is.na(header_idx)) return(NA_real_)
  header_fields <- strsplit(lines[header_idx], "\t", fixed = TRUE)[[1]]
  data_lines <- lines[(header_idx + 1):length(lines)]
  
# Split each data line on tab
  split_data <- strsplit(data_lines, "\t", fixed = TRUE)
  max_len <- max(lengths(split_data))
  
  mat <- t(vapply(
    split_data,
    function(x) c(x, rep(NA_character_, max_len - length(x))),
    character(max_len) ))
  
  dat <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(dat) <- header_fields
  
#need the gene_name column to find NKX2-1
  if (!"gene_name" %in% colnames(dat)) return(NA_real_)
  nkx_row <- dat[dat$gene_name == "NKX2-1", , drop = FALSE]
  if (nrow(nkx_row) == 0) return(NA_real_)
  as.numeric(nkx_row$tpm_unstranded[1])
}


## 4. Apply function to all samples and compute log2(TPM + 1)
ss_rna$TPM_NKX2_1 <- sapply(ss_rna$Dir, get_tpm_NKX2)

ss_clean <- subset(ss_rna, !is.na(TPM_NKX2_1))
ss_clean$log2_TPM_plus_1 <- log2(ss_clean$TPM_NKX2_1 + 1)


## 5. Plot tumor vs normal NKX2-1 expression
boxplot(
  log2_TPM_plus_1 ~ Sample.Type,
  data = ss_clean,
  main = "NKX2-1 Expression in LUAD (TCGA)",
  xlab = "Tissue Type",
  ylab = "log2(TPM + 1)",
  col  = c("tomato", "skyblue")
)
