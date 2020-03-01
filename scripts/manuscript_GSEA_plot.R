##############################
# manuscript_GSEA_plot.R
##############################

# Show GSEA results of GSEA diff expression

# Load GSEA results 
cond_names<- c("Ctrl+RA Vs. Ctrl","DLG2 Vs. Ctrl","DLG2+RA Vs. Ctrl")
GSEA_ls<- readRDS("results/data/GSEA_results.rds")
for(i in 1:2){
  if(i==1) pdf(paste0("results/figs/manuscript_GSEA.pdf"))
  else svglite::svglite("results/figs/manuscript_GSEA.svg")
  par(mfrow=c(3,1))
  for(j in 1:3){
    GSEA<- GSEA_ls[[j]][[1]] # [[1]] = GO   
    if(j==1) cond_temp<- "RA"
    if(j==2) cond_temp<- "DLG2"
    if(j==3) cond_temp<- "DLG2_RA"
    
    # Visualize top 5 diff related GO
    GSEA_diff<- GSEA[grepl("DIFFERENTIATION",rownames(GSEA)),] # 196
    barplot(-log10(as.numeric(as.character(GSEA_diff[5:1,"q"]))),main=cond_temp,horiz = T)
    text(rep(0,10),seq(11.5,0.7,by = -1.2),labels = rownames(GSEA_diff)[1:5],adj = 0)  
  }
  dev.off()  
}

