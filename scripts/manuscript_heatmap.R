##############################
# manuscript_heatmap.R
##############################

# Aim: show heatmap comparing diff conditions

library(gplots) # Heatmap 2

# Load diff expressed genes
results_diff_expr<- readRDS("data/DLG2_diff_expression_results.rds")[1:3] # Only diff expr to control

# Only focus on selected genes: 1071 in total (798 with expr filter)
selected_genes<- NULL
for(comparison in names(results_diff_expr)){
  res_proc<- results_diff_expr[[comparison]]
  res_proc<- as.data.frame(res_proc)
  res_proc<- res_proc[res_proc$baseMean>10,] # Select minimal expression
  idx_sign<- which(-log10(res_proc[,"padj"])>2) 
  idx_up<- which(res_proc[,"log2FoldChange"]>2)
  idx_down<- which(res_proc[,"log2FoldChange"]<(-2)) 
  idx_diff_sign<- intersect(idx_sign,union(idx_up,idx_down)) 
  genes_sel<- rownames(res_proc[idx_diff_sign,])
  selected_genes<- unique(c(selected_genes,genes_sel))
}

# Get matrix of selected genes
logFC_matrix<- NULL
for(comparison in names(results_diff_expr)){
  res_proc<- results_diff_expr[[comparison]]
  res_proc<- as.data.frame(res_proc)
  logFC_matrix<- cbind(logFC_matrix,res_proc[selected_genes,"log2FoldChange"])
}
rownames(logFC_matrix)<- selected_genes
colnames(logFC_matrix)<- names(results_diff_expr)

# Set cut-offs for visualization
data<- logFC_matrix
data[data>=3]<-3
data[data<=-3]<--3

# Plot dummy heatmap to define clusters afterwars
pdf(file = NULL)
hm<- heatmap.2(data)
# hm<- heatmap.2(data,hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off()

# Get main clusters: set to 8 to define but redefine later
n_classes<- 4
clust_classes<- cutree(as.hclust(hm$rowDendrogram),n_classes)

# Get matrix coordinates of sepline
sepline<- cumsum(table(factor(clust_classes[hm$rowInd],levels = rev(unique(clust_classes[hm$rowInd])))))

for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_heatmap.pdf")
  else svglite::svglite("results/figs/manuscript_heatmap.svg")
  # hm<- heatmap.2(data,rowsep =sepline,RowSideColors= rainbow(n_classes)[clust_classes],sepcolor="black",sepwidth=c(3,3),col=bluered(75),symkey=TRUE, key=TRUE, keysize=1, trace="none",scale="none",density.info='none', cexRow=0.7, cexCol=0.5,key.title = "logFC")
  hm<- heatmap.2(data[,c(2,3,1)],rowsep =sepline,RowSideColors= rainbow(n_classes)[clust_classes],sepcolor="black",sepwidth=c(3,3),col=bluered(75),symkey=TRUE, key=TRUE, keysize=1, trace="none",scale="none",density.info='none', cexRow=0.7, cexCol=0.5,key.title = "logFC",Colv = F,dendrogram="row")
  # hm<- heatmap.2(data,rowsep =sepline,RowSideColors= rainbow(n_classes)[clust_classes],sepcolor="black",sepwidth=c(3,3),col=bluered(75),symkey=TRUE, key=TRUE, keysize=1, trace="none",scale="none",density.info='none', cexRow=0.7, cexCol=0.5,key.title = "logFC",hclustfun=function(x) hclust(x, method="ward.D2"))
  legend("topright",legend = paste("Cluster", 1:n_classes),fill = rainbow(n_classes),cex=0.5)
  dev.off()
}

# GSEA TF
source("scripts/functions/do_GSEA2.R")
# load("downloads/TFT/TFT.RData")
TFT_RN<- readRDS("downloads/regNetwork/TFT_RN.rds")
genes_all<- rownames(results_diff_expr[[1]])
TF_sign<- NULL
# TFT_RN<- regulon_ls
for(i in 1:n_classes){
  genes_sel<- names(rev(clust_classes[hm$rowInd])[rev(clust_classes[hm$rowInd])==i])
  GSEA_RN<- do_GSEA2(genes_retrieved = genes_sel,genes_all = genes_all,GSEA_db = TFT_RN,min_genes = 4,isList = T)
  # GSEA_RN<- do_GSEA2(genes_retrieved = genes_sel,genes_all = genes_all,GSEA_db = TFT_TRRUST2,min_genes = 4,isList = T)
  assign(paste0("GSEA_RN_Clust_",i),GSEA_RN)  
  TF_sign_temp<- rownames(GSEA_RN)[as.numeric(as.character(GSEA_RN[,"q"]))<0.1]
  TF_sign<- rev(sort(unique(c(TF_sign,TF_sign_temp))))
  # cat("Cluster",i,":",rownames(GSEA_RN)[as.numeric(as.character(GSEA_RN[,"q"]))<0.1],"\n")
}

# Cluster 1 : RARA RARB RELA AR 
# Cluster 2 : STAT5B CEBPB 
# Cluster 3 : POU2F1 MEF2A 
# Cluster 4 : REST TBP SMAD1 TCF4 

# All expressed? no AR! exclude
DLG2_normalized_counts <- readRDS("data/DLG2_normalized_counts.rds")
DLG2_normalized_counts[TF_sign,]
TF_sign<- TF_sign[TF_sign!="AR"]
  
for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_heatmap_TF.pdf")
  else svglite::svglite("results/figs/manuscript_heatmap_TF.svg")
  par(mfrow=c(2,4))
  clust_TF_matrix<- matrix(NA,length(TF_sign),n_classes,dimnames = list(TF_sign,1:n_classes))
  for(i in 1:n_classes){
    GSEA_RN<- get(paste0("GSEA_RN_Clust_",i))
    # barplot(-log10(as.numeric(as.character(GSEA_RN[rev(TF_sign),"q"]))),main=i,horiz = T,names.arg = rev(TF_sign),las=2,xlim=c(0,10))
    bg_prop<- 100*as.numeric(as.character(GSEA_RN[rev(TF_sign),"n_genes_pw"]))/length(genes_all)
    pw_prop<- as.numeric(as.character(GSEA_RN[rev(TF_sign),"prop_pos"]))   
    isSign<- as.numeric(as.character(GSEA_RN[rev(TF_sign),"q"]))<0.1   
    barplot(pw_prop/bg_prop,main=i,horiz = T,names.arg = rev(TF_sign),las=2,xlim=c(0,12),xlab="Times enriched compared to background",col=isSign)
    abline(v=1,lty=2)
    clust_TF_matrix[rev(TF_sign),i]<- pw_prop/bg_prop
  }
  dev.off()
}

# 1 plot
for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_heatmap_TF_1plot.pdf")
  else svglite::svglite("results/figs/manuscript_heatmap_TF_1plot.svg")
  par(mfrow=c(1,4))
  barplot(t(clust_TF_matrix[,c(3,4,1,2)]),beside=T,horiz = T,las=2,col=rainbow(n_classes)[c(3,4,1,2)],xlab="Times enriched")
  abline(v=1,lty=2)
  dev.off()
}
