##########################
# manuscript_bridge.R
##########################

# Aim: correlate Bridge genes (Furlan et al.) to NG prognosis (R2 platform)

# Load & process data
#####################
{
  # Load Bridge genes
  # See fig 6A from manuscript
  # Genes given in table S1 (E12.5), although many more seem to be present
  table_S1 <- readxl::read_excel("downloads/pub/furlan_2017/NIHMS974217-supplement-Tables_S1-S2.xlsx")[,1:4]
  # Manually check
  bridge_genes<-c(
    'Igfbp3',
    'Lrig1',
    'Tubb6',
    'Cyp2j9',
    'E130114P18Rik',
    'Samd5',
    'Nrarp',
    'Fndc3c1',
    'Nme4',
    'Sox2',
    'Hmcn1',
    'Slc39a8',
    'Rgs16',
    'Tpx2',
    'Igfbp2',
    'Shmt1',
    'Chek1',
    'Prim1',
    'Arvcf',
    'Bub1',
    'Cdca8',
    'Ncapg',
    'Spag5',
    'Tyms',
    'Thbd',
    'Kit',
    'Ier5',
    'Tbx2',
    'Shcbp1',
    'Ckb',
    'Mns1',
    'Pbk',
    'Cenpe',
    'Kntc1',
    'Jam2',
    'Pde5a',
    'Tk1',
    'Top2a',
    'Id2',
    'Cenpf',
    'Csrp2',
    'Spc25',
    'Ascl1',
    'Cdkn1c',
    'Rftn2',
    'Fanca',
    'Tlx2',
    'Fat1',
    'Cenpk',
    'Dqx1',
    'Rrm2',
    'Ube2c',
    'Cenpm',
    'Sct',
    'Phox2b',
    'Melk',
    'Sox11',
    'Htr3b',
    'Htr3a',
    'Rtn4rl1',
    'Dll3',
    'Mfng',
    'Gadd45g',
    'Igsf9',
    'Gcnt2',
    'Hes6',
    'Plxna2',
    'Fam110a',
    'Setbp1',
    'Phf21b',
    'Pdlim3',
    'Tbx20',
    'Cttnbp2',
    'Nefl',
    'Nefm',
    'Fbrsl1',
    'Kcnj12',
    'Tagln3',
    'Stra6',
    'Hipk2',
    'Brsk2',
    'Efnb2',
    'Slc10a4',
    'Cdc25b',
    'Pde2a',
    'Dpysl3',
    'Grik3',
    'Mapk11',
    'Ank',
    'Miat',
    'Glt28d2',
    # 'Adcy', # Not sure whether they mean or Adcyap1r1 or Adcy1? Exclude
    'Msrb2',
    'C1ql1',
    'Dlgap2',
    'Igsf21',
    'Sept3',
    'Adamts19',
    'Dlg2',
    'Tmem179'
  )
  
  # Load R2 platform genes
  # Data not downloadable? 
  # Manual curation performed by Bengt & Joakim + double check afterwards by Joachim
  R2<- cbind(bridge_genes,c(
    0,
    1,
    1,
    NA,
    NA,
    1,
    NA,
    NA,
    0,
    0,
    1,
    0,
    1,
    0,
    0,
    0,
    0,
    0,
    1,
    0,
    0,
    0,
    NA,
    0,
    NA,
    NA,
    NA,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1,
    1,
    0,
    0,
    1,
    0,
    1,
    0,
    0,
    0,
    0,
    1,
    0,
    0,
    0,
    1,
    0,
    0,
    0,
    1,
    0,
    1,
    0,
    0,
    0,
    NA,
    0,
    0,
    1,
    NA,
    1,
    1,
    1,
    NA,
    NA,
    0,
    0,
    1,
    1,
    0,
    0,
    1,
    1,
    1,
    1,
    NA,
    1,
    NA,
    # 1,
    1,
    0,
    1,
    1,
    1,
    0,
    1,
    1
  ))
  
  # Add HGNC names
  library(biomaRt)
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl") 
  map_MGI_HGNC<- getLDS(attributes = c("hgnc_symbol"), 
                        filters = "hgnc_symbol", values = bridge_genes, mart = human, 
                        attributesL = c("mgi_symbol"), martL = mouse)
  map_MGI_HGNC<- map_MGI_HGNC[map_MGI_HGNC$MGI.symbol%in%bridge_genes,]
  # Add some manually
  unmapped<- setdiff(bridge_genes,map_MGI_HGNC$MGI.symbol)
  # "Cyp2j9"        "E130114P18Rik" "Fndc3c1"       "Nefm"         "Ank"           "Miat"          "Glt28d2"      "Sept3" 
  map_MGI_HGNC<- rbind(map_MGI_HGNC,cbind(HGNC.symbol=c(NA,NA,NA,"NEFM","ANKH","MIAT",NA,"SEPT3"),MGI.symbol=unmapped))

  # Peak time & expr from table S2
  peak_time<- table_S1$peak.time
  names(peak_time)<- table_S1$gene
  peak_expr<- table_S1$peakM
  names(peak_expr)<- table_S1$gene
  
  # Fuse datasets
  rownames(R2)<- R2[,1]
  bridge_matrix<- cbind(map_MGI_HGNC,R2[map_MGI_HGNC$MGI.symbol,2],peak_time[map_MGI_HGNC$MGI.symbol],peak_expr[map_MGI_HGNC$MGI.symbol],stringsAsFactors=F)
  colnames(bridge_matrix)<- c("HGNC","MGI","hasGoodProgn","peak_time","peak_expr")
  bridge_matrix<- na.omit(bridge_matrix)
  bridge_matrix<- bridge_matrix[order(bridge_matrix$peak_time),]
  bridge_matrix$hasGoodProgn<- as.numeric(bridge_matrix$hasGoodProgn)
}

# 1) Correlation with prognosis
################################

# Logistic regression on peak time
svglite::svglite("results/figs/bridge_peak Vs_R2.svg")
glm_model<- glm(hasGoodProgn~peak_time,data = bridge_matrix,family=binomial(link = "logit"))
pred1<-predict(glm_model,data.frame(peak_time=bridge_matrix$peak_time),type="response",se.fit=TRUE)
plot(bridge_matrix$peak_time,100*pred1$fit,ylim=c(0,100),type="l",frame.plot=F,xlab="Peak time",ylab="% Good prognosis (R2)",col="blue",lwd=3)
polygon(c(bridge_matrix$peak_time,rev(bridge_matrix$peak_time)), c(100*pred1$fit+196*pred1$se.fit, rev(100*pred1$fit-196*pred1$se.fit)),col = rgb(.9,.9,.9,.5), border = NA)
points(bridge_matrix$peak_time,100*as.numeric(bridge_matrix$hasGoodProgn),pch=16)
glm_model$coefficients["peak_time"] # 4.85
summary(glm_model)$coefficients["peak_time",4] # p=0.001920145
# Label
labelled_genes<- c("ASCL1", "PHOX2B", "CHEK1", "NEFM", "SOX2", "IGFBP2", "ID2", "DLGAP2", "IGSF21", "ADAMTS19", "DLG2", "TMEM179")
bridge_matrix_labelled<- bridge_matrix[bridge_matrix$HGNC%in%labelled_genes,]
points(bridge_matrix_labelled$peak_time,100*as.numeric(bridge_matrix_labelled$hasGoodProgn),pch=16,col="blue")
text(bridge_matrix_labelled$peak_time,100*as.numeric(bridge_matrix_labelled$hasGoodProgn),labels = bridge_matrix_labelled$HGNC,srt=90,adj=c(0,0.5),xpd=NA,cex=0.7)
dev.off()

# Save data
#############
save(bridge_matrix, map_MGI_HGNC,file="results/data/bridge_analysis.RData")



