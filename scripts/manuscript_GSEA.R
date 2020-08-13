##############################
# manuscript_GSEA.R
##############################

source("scripts/functions/do_GSEA2.R")
load("downloads/mSigDB/v62/MSigDB_v62.RData")
results_diff_expr<- readRDS("data/DLG2_diff_expression_results.rds")

GSEA_results_ls<- list()
for(comparison in names(results_diff_expr)){
  res_proc<- results_diff_expr[[comparison]]
  res_proc<- as.data.frame(res_proc)
  idx_sign<- which(-log10(res_proc[,"padj"])>2) # 7303 significant!!!
  idx_up<- which(res_proc[,"log2FoldChange"]>2) # 1123
  idx_down<- which(res_proc[,"log2FoldChange"]<(-2)) # 384
  idx_diff_sign<- intersect(idx_sign,union(idx_up,idx_down)) # 650
  genes_sel<- rownames(res_proc[idx_diff_sign,])
  genes_all<- rownames(res_proc)
  
  # GO: lot ECM, also migration, ...
  GSEA_GO<- do_GSEA2(genes_retrieved = genes_sel,genes_all = genes_all,GSEA_db = GO_ls,min_genes = 2,isList = T)
  
  # Rea: GPCR, rhodopsin, ...
  GSEA_Rea<- do_GSEA2(genes_retrieved = genes_sel,genes_all = genes_all,GSEA_db = Rea_ls,min_genes = 2,isList = T)
  
  # Kegg: neuroactive ligands, ECM, ...
  GSEA_Kegg<- do_GSEA2(genes_retrieved = genes_sel,genes_all = genes_all,GSEA_db = Kegg_ls,min_genes = 2,isList = T)
  
  # Put results in list
  GSEA_results_ls<- c(GSEA_results_ls,list(list(GSEA_GO,GSEA_Rea,GSEA_Kegg)))
}
names(GSEA_results_ls)<- names(results_diff_expr)

# Save
saveRDS(GSEA_results_ls,file = "results/data/GSEA_results.rds")

