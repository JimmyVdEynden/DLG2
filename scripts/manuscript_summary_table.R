# manuscript_summary_table.R
##########################

# Aim: summarize results in 1 suppl. table

# RNA-Seq
DLG2_diff_expression_results <- readRDS("data/DLG2_diff_expression_results.rds")
GSEA_results <- readRDS("results/data/GSEA_results.rds")

# Proteomics
DLG2_diff_expression_proteomics_results <- readRDS("data/DLG2_diff_expression_proteomics_results.rds")

# Diff markers
load("downloads/pub/van_groningen_2017/vGron_diff_markers.RData")

# Fuse in 1 table
for(i in c(2,3,1)){
  # RNA-Seq
  DE_temp<- as.data.frame(DLG2_diff_expression_results[[i]])
  DE_temp<- DE_temp[,c("log2FoldChange","pvalue")]
  colnames(DE_temp)<- c("RNA_logFC","RNA_p")
  # Prot
  DE_prot_temp<- DLG2_diff_expression_proteomics_results[[i]]
  DE_prot_temp<- as.data.frame(DE_prot_temp[,c("logFC","pvalue")])
  colnames(DE_prot_temp)<- c("PROT_logFC","PROT_p")
  # Combine
  common_genes<- sort(union(rownames(DE_temp),rownames(DE_prot_temp)))
  DE_comb<- cbind(DE_temp[common_genes,],DE_prot_temp[common_genes,])
  # Add diff markers
  DE_comb<- cbind(DE_comb,diff_marker=NA)
  DE_comb$diff_marker[rownames(DE_comb)%in%ADRN]<- "ADRN"
  DE_comb$diff_marker[rownames(DE_comb)%in%MES]<- "MES"
  # Assign
  if(i==1) DE_comb_RA<- DE_comb  
  if(i==2) DE_comb_DLG2<- DE_comb  
  if(i==3) DE_comb_DLG2_RA<- DE_comb  
}

# Save to excel
WriteXLS::WriteXLS(c("DE_comb_DLG2","DE_comb_DLG2_RA","DE_comb_RA"), "results/tables/manuscript_summary.xlsx",SheetNames = c("DLG2","RA+DLG2","RA"),row.names = T)
