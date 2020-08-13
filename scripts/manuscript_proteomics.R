#################################
# manuscript_proteomics.R
#################################

# Show diff expression vulcano plot of DLG2 wity expression markers

# Load diff expression results
results_diff_expr<- readRDS(file = "data/DLG2_diff_expression_proteomics_results.rds")

# Get DLG2 results
res_proc<- results_diff_expr[["DLG2 Vs. Ctrl"]]

# Plot Volcano
###############
for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_vulcano_proteomics.pdf")
  else svglite::svglite("results/figs/manuscript_vulcano_proteomics.svg")
  # res_proc<- res_proc[!is.na(res_proc[,"pvalue"]),]
  idx_sign<- which(-log10(res_proc[,"pvalue"])>2) # 691 significant
  idx_up<- which(res_proc[,"logFC"]>1) # 87
  idx_down<- which(res_proc[,"logFC"]<(-1)) # 253
  idx_up_sign<- intersect(idx_sign,idx_up) # 49
  idx_down_sign<- intersect(idx_sign,idx_down) # 30
  
  # par(mfrow=c(1,2))
  plot(res_proc[-c(idx_up_sign,idx_down_sign),"logFC"],-log10(res_proc[-c(idx_up_sign,idx_down_sign),"pvalue"]),pch=16,col=rgb(0.8,0.8,0.8,0.5),frame.plot=F,xlab="log2FC",ylab="-log10(Pvalue)",xlim=c(-4,4),xpd=NA,main="DLG2 overexpression")
  abline(h=2,v=c(-1,1),lty=2)
  points(res_proc[idx_up_sign,"logFC"],-log10(res_proc[idx_up_sign,"pvalue"]),pch=16,col=rgb(0,0,1,0.5))
  points(res_proc[idx_down_sign,"logFC"],-log10(res_proc[idx_down_sign,"pvalue"]),pch=16,col=rgb(0,0,1,0.5))
  points(res_proc["DLG2","logFC"],-log10(res_proc["DLG2","pvalue"]),cex=2,pch=16,col="purple")
  text(res_proc["DLG2","logFC"],-log10(res_proc["DLG2","pvalue"]),"DLG2",adj = c(0,0),xpd=NA,col="purple")  
  
  # Plot diff markers
  diff_markers<- c("NEFM","NEFL","NEFH","RET","NRTK1","MAP2","RBFOX3","NRCAM","ASCL1","MYCN","VGF", "SGK1", "EGR1", "NTRK1", "DLGAP2") # Manually defined
  idx_diff_sign_marker<- which(rownames(res_proc)%in%diff_markers)
  points(res_proc[idx_diff_sign_marker,"logFC"],-log10(res_proc[idx_diff_sign_marker,"pvalue"]),pch=16,col=rgb(1,0,0,1))
  for(gene in diff_markers){
    if(!gene%in%rownames(res_proc)) next
    text(res_proc[gene,"logFC"],-log10(res_proc[gene,"pvalue"]),gene,adj = c(0,0),xpd=NA,col="red")  
  }
  legend("topleft",legend = c(paste0(length(idx_up_sign)," upregulated"),paste0(length(idx_down_sign))," downregulated"),bty="n",text.col="blue")
  dev.off()  
}

# heatmap
##########

library(gplots) # Heatmap 2

# Only focus on selected genes: 117
selected_genes<- NULL
for(comparison in names(results_diff_expr)){
  res_proc<- results_diff_expr[[comparison]]
  idx_sign<- which(-log10(res_proc[,"pvalue"])>2) 
  idx_up<- which(res_proc[,"logFC"]>1)
  idx_down<- which(res_proc[,"logFC"]<(-1)) 
  idx_diff_sign<- intersect(idx_sign,union(idx_up,idx_down)) 
  genes_sel<- rownames(res_proc[idx_diff_sign,])
  selected_genes<- unique(c(selected_genes,genes_sel))
}

# Get matrix of selected genes
logFC_matrix<- NULL
for(comparison in names(results_diff_expr)){
  res_proc<- results_diff_expr[[comparison]]
  logFC_matrix<- cbind(logFC_matrix,res_proc[selected_genes,"logFC"])
}
rownames(logFC_matrix)<- selected_genes
colnames(logFC_matrix)<- names(results_diff_expr)

# Set cut-offs for visualization
data<- logFC_matrix[,c(2,3,1)]
data[data>=3]<-3
data[data<=-3]<--3

# Plot dummy heatmap to define clusters afterwars
pdf(file = NULL)
hm<- heatmap.2(data)
dev.off()

# Get main clusters:
n_classes<- 4
clust_classes<- cutree(as.hclust(hm$rowDendrogram),n_classes)

# Get matrix coordinates of sepline
sepline<- cumsum(table(factor(clust_classes[hm$rowInd],levels = rev(unique(clust_classes[hm$rowInd])))))

# Plot
for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_heatmap_proteomics.pdf")
  else svglite::svglite("results/figs/manuscript_heatmap_proteomics.svg")
  hm<- heatmap.2(data,rowsep =sepline,RowSideColors= rainbow(n_classes)[clust_classes],sepcolor="black",sepwidth=c(0.5,0.5),col=bluered(75),symkey=TRUE, key=TRUE, keysize=1, trace="none",scale="none",density.info='none', cexRow=0.5, cexCol=0.5,key.title = "logFC",Colv = F,dendrogram="row")
  legend("topright",legend = paste("Cluster", 1:n_classes),fill = rainbow(n_classes),cex=0.5)
  dev.off()
}
