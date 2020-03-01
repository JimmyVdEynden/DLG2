#################################
# manuscript_vulcano.R
#################################

library(gplots)

# Show diff expression vulcano plot wity differentiation markers

# Load diff expression results
results_diff_expr<- readRDS(file = "data/DLG2_diff_expression_results.rds")[1:3]
names(results_diff_expr)<- c("RA","DLG2","DLG2+RA")

# Differentiation markers to label
diff_markers<- c("NEFM","NEFL","NEFH","RET","NRTK1","MAP2","RBFOX3","NRCAM","ASCL1","MYCN","VGF", "SGK1", "EGR1", "NTRK1", "DLGAP2") # Manually defined

for(cond in names(results_diff_expr)){
  
  # Get DE results
  res_proc<- results_diff_expr[[cond]]
  res_proc<- as.data.frame(res_proc)
  
  # Plot
  for(i in 1:2){
    if(i==1) pdf(paste0("results/figs/manuscript_vulcano_",cond,".pdf"))
    else svglite::svglite(paste0("results/figs/manuscript_vulcano_",cond,".svg"))
    res_proc<- res_proc[!is.na(res_proc[,"padj"]),] # remove NA (i.e. no expression before/after treatment), 17242 left
    res_proc[res_proc[,"padj"]==0,"padj"]<- 1e-300 # if q=0, put at lowest possible value
    idx_sign<- which(-log10(res_proc[,"padj"])>2) # 5495 significant
    idx_up<- which(res_proc[,"log2FoldChange"]>2) # 917
    idx_down<- which(res_proc[,"log2FoldChange"]<(-2)) # 253
    idx_up_sign<- intersect(idx_sign,idx_up) # 414
    idx_down_sign<- intersect(idx_sign,idx_down) # 122
    idx_diff_sign<- c(idx_up_sign,idx_down_sign)
    
    # Genes
    genes_DE<- rownames(res_proc)[c(idx_down_sign,idx_up_sign)]
    assign(paste0("genes_",cond),genes_DE)
    
    # par(mfrow=c(1,2))
    plot(res_proc[-idx_diff_sign,"log2FoldChange"],-log10(res_proc[-idx_diff_sign,"padj"]),pch=16,col=rgb(0.8,0.8,0.8,0.5),frame.plot=F,xlab="log2FC",ylab="-log10(FDR)",xlim=c(-12,12),ylim=c(0,300),xpd=NA,main=cond)
    abline(h=2,v=c(-2,2),lty=2)
    text(12,2,"-log10(FDR)=2",adj=c(1,0))
    text(-2,300,"log2FC=-2",adj=c(1,0),srt=90)
    text(2,300,"log2FC=2",adj=c(1,0),srt=90)
    points(res_proc[idx_up_sign,"log2FoldChange"],-log10(res_proc[idx_up_sign,"padj"]),pch=16,col=rgb(0,0,1,0.5))
    points(res_proc[idx_down_sign,"log2FoldChange"],-log10(res_proc[idx_down_sign,"padj"]),pch=16,col=rgb(0,0,1,0.5))
    points(res_proc["DLG2","log2FoldChange"],-log10(res_proc["DLG2","padj"]),cex=2,pch=16,col="purple")
    text(res_proc["DLG2","log2FoldChange"],-log10(res_proc["DLG2","padj"]),"DLG2",adj = c(0,0),xpd=NA,col="purple")  
    
    # Plot diff markers
    idx_diff_sign_marker<- which(rownames(res_proc)%in%diff_markers)
    points(res_proc[idx_diff_sign_marker,"log2FoldChange"],-log10(res_proc[idx_diff_sign_marker,"padj"]),pch=16,col=rgb(1,0,0,1))
    for(gene in diff_markers){
      if(!gene%in%rownames(res_proc)) next
      text(res_proc[gene,"log2FoldChange"],-log10(res_proc[gene,"padj"]),gene,adj = c(0,0),xpd=NA,col="red")  
    }
    
    # Legend
    legend("topleft",legend = c(paste0(length(idx_up_sign)," upregulated"),paste0(length(idx_down_sign))," downregulated"),bty="n",text.col="blue")
    dev.off()  
  }
  
}

# Venn diagram
##############
venn(list(genes_DLG2,genes_RA,`genes_DLG2+RA`))
venn(list(genes_DLG2,genes_RA))

# GSEA TF
source("scripts/functions/do_GSEA2.R")
# load("downloads/TFT/TFT.RData")
TFT_RN<- readRDS("downloads/regNetwork/TFT_RN.rds")
genes_all<- rownames(results_diff_expr[[1]])
genes_sel<- intersect(intersect(genes_RA,`genes_DLG2+RA`),genes_DLG2)

# Only expressed TFs
DLG2_normalized_counts <- readRDS("data/DLG2_normalized_counts.rds")
genes_expr<- rownames(DLG2_normalized_counts)[rowSums(DLG2_normalized_counts>=10)>0]
TFT_RN_expr<- TFT_RN[names(TFT_RN)%in%genes_expr]

GSEA_RN<- do_GSEA2(genes_retrieved = genes_sel,genes_all = genes_all,GSEA_db = TFT_RN_expr,min_genes = 4,isList = T)






