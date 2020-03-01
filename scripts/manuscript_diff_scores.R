##############################
# manuscript_diff_scores.R
##############################

# Aim: show MES and ADRN differentiation scores for different conditions

# Get expression counts
DLG2_normalized_counts <- readRDS("data/DLG2_normalized_counts.rds")

# Get diff expr results
results_diff_expr<- readRDS(file = "data/DLG2_diff_expression_results.rds")

# Load MES and ADRN diff markers as published by V. Groningen et al.
# load("data/diff_markers.RData")
load("downloads/pub/van_groningen_2017/vGron_diff_markers.RData")

# 1) Determine enrichment marker
#################################

# ADRN
ADRN_enrich_res<- matrix(NA,3,2,dimnames = list(c("RA","DLG2","DLG2+RA"),c("enrich","p")))
for(i in 1:nrow(ADRN_enrich_res)){
  res_proc<- results_diff_expr[[i]]
  res_proc<- as.data.frame(res_proc)
  all_genes<- rownames(res_proc)
  idx_sign<- which(-log10(res_proc[,"padj"])>2) 
  idx_diff<- which(abs(res_proc[,"log2FoldChange"])>2)
  idx_diff_sign<- intersect(idx_sign,idx_diff) 
  isADRN<- all_genes%in%ADRN
  isDE<- all_genes%in%rownames(res_proc)[idx_diff_sign]
  ADRN_DE_t<- table(isADRN,isDE)
  ADRN_DE_prop<- prop.table(ADRN_DE_t,2)[2,]
  ADRN_DE_enrich<- ADRN_DE_prop["TRUE"]/ADRN_DE_prop["FALSE"]
  ADRN_DE_t_ft<- fisher.test(ADRN_DE_t,alternative = "greater")
  ADRN_enrich_res[i,]<- c(ADRN_DE_enrich,ADRN_DE_t_ft$p.value)
}

# MES
MES_enrich_res<- matrix(NA,3,2,dimnames = list(c("RA","DLG2","DLG2+RA"),c("enrich","p")))
for(i in 1:nrow(MES_enrich_res)){
  res_proc<- results_diff_expr[[i]]
  res_proc<- as.data.frame(res_proc)
  all_genes<- rownames(res_proc)
  idx_sign<- which(-log10(res_proc[,"padj"])>2) 
  idx_diff<- which(abs(res_proc[,"log2FoldChange"])>2)
  idx_diff_sign<- intersect(idx_sign,idx_diff) 
  isMES<- all_genes%in%MES
  isDE<- all_genes%in%rownames(res_proc)[idx_diff_sign]
  MES_DE_t<- table(isMES,isDE)
  MES_DE_prop<- prop.table(MES_DE_t,2)[2,]
  MES_DE_enrich<- MES_DE_prop["TRUE"]/MES_DE_prop["FALSE"]
  MES_DE_t_ft<- fisher.test(MES_DE_t,alternative = "greater")
  MES_enrich_res[i,]<- c(MES_DE_enrich,MES_DE_t_ft$p.value)
}

# Plot
for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_diff_marker_enrichment.pdf")
  else svglite::svglite("results/figs/manuscript_diff_marker_enrichment.svg")
  par(mfrow=c(1,3))
  # ADRN
  bp<- barplot(ADRN_enrich_res[c("DLG2","DLG2+RA","RA"),1],ylab="Marker gene enrichment",ylim=c(0,3),main="ADRN")
  abline(h=1,lty=2)
  text(bp,c(3,3,3),paste0("P=",format(ADRN_enrich_res[c("DLG2","DLG2+RA","RA"),2],scientific = T,digits = 3)),xpd=NA)
  # MES
  bp<- barplot(MES_enrich_res[c("DLG2","DLG2+RA","RA"),1],ylab="Marker gene enrichment",ylim=c(0,3),main="MES")
  abline(h=1,lty=2)
  text(bp,c(3,3,3),paste0("P=",format(MES_enrich_res[c("DLG2","DLG2+RA","RA"),2],scientific = T,digits = 3)),xpd=NA)
  dev.off()
}

# 2) Calulate scores
#####################

# Get Rank score for all genes
DLG2_scores<- apply(DLG2_normalized_counts, 2, function(x) rank(x)/length(x))

# Get MES and ADRN scores
MES_scores<- DLG2_scores[rownames(DLG2_scores)%in%MES,]
ADRN_scores<- DLG2_scores[rownames(DLG2_scores)%in%ADRN,]

# Compare diff scores cell lines
for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_diff_scores.pdf")
  else svglite::svglite("results/figs/manuscript_diff_scores.svg")
  plot(colMeans(ADRN_scores),colMeans(MES_scores),pch=16,col=rep(c("blue","purple","yellow","red"),each=4),xlab="ADRN Score",ylab="MES Score",cex=2)
  legend("topright",legend = c("Ctrl","Ctrl_RA","DLG2","DLG2_RA"),col=rep(c("yellow","red","blue","purple")),pch=16)
  # Add mean
  DLG2_scores_mean<- cbind(rowMeans(DLG2_scores[,1:4]),rowMeans(DLG2_scores[,5:8]),rowMeans(DLG2_scores[,9:12]),rowMeans(DLG2_scores[,13:16]))
  colnames(DLG2_scores_mean)<- c("DLG2","DLG2_RA","Ctrl","Ctrl_RA")
  rownames(DLG2_scores_mean)<- rownames(DLG2_normalized_counts)
  MES_scores_mean<- DLG2_scores_mean[rownames(DLG2_scores_mean)%in%MES,]
  ADRN_scores_mean<- DLG2_scores_mean[rownames(DLG2_scores_mean)%in%ADRN,]
  points(colMeans(ADRN_scores_mean),colMeans(MES_scores_mean),pch="+",col=c("blue","purple","yellow","red"),cex=2)
  dev.off()
}


# Significance? Decide later how to plot
aov_ADRN<- aov(colMeans(ADRN_scores)~factor(rep(c("DLG2","DLG2_RA","Ctrl","Ctrl_RA"),each=4)))
summary(aov_ADRN)
TukeyHSD(aov_ADRN)
# Ctrl_RA-Ctrl    -0.017423174 -0.020897128 -0.01394922 0.000000
# DLG2-Ctrl        0.026130110  0.022656156  0.02960406 0.000000
# DLG2_RA-Ctrl     0.004963736  0.001489781  0.00843769 0.005437
# DLG2-Ctrl_RA     0.043553284  0.040079330  0.04702724 0.000000
# DLG2_RA-Ctrl_RA  0.022386910  0.018912955  0.02586086 0.000000
# DLG2_RA-DLG2    -0.021166375 -0.024640329 -0.01769242 0.000000
t.test(colMeans(ADRN_scores)[9:12],colMeans(ADRN_scores)[1:4]) # p-value = 1.232e-06
t.test(colMeans(ADRN_scores)[9:12],colMeans(ADRN_scores)[5:8]) # p-value = 0.001067
t.test(colMeans(ADRN_scores)[9:12],colMeans(ADRN_scores)[13:16]) # p-value = 0.0001456
aov_MES<- aov(colMeans(MES_scores)~factor(rep(c("DLG2","DLG2_RA","Ctrl","Ctrl_RA"),each=4)))
summary(aov_MES)
TukeyHSD(aov_MES)
# Ctrl_RA-Ctrl     0.024549582  0.016629966  0.032469198 0.0000046
# DLG2-Ctrl       -0.002284344 -0.010203960  0.005635272 0.8267209
# DLG2_RA-Ctrl     0.026568404  0.018648788  0.034488020 0.0000020
# DLG2-Ctrl_RA    -0.026833926 -0.034753542 -0.018914310 0.0000018
# DLG2_RA-Ctrl_RA  0.002018822 -0.005900794  0.009938438 0.8720635
# DLG2_RA-DLG2     0.028852748  0.020933132  0.036772364 0.0000008
t.test(colMeans(MES_scores)[9:12],colMeans(MES_scores)[1:4]) # p-value = 0.3234
t.test(colMeans(MES_scores)[9:12],colMeans(MES_scores)[5:8]) # p-value = 1.857e-05
t.test(colMeans(MES_scores)[9:12],colMeans(MES_scores)[13:16]) # p-value = 0.001678


