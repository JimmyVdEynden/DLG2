# Download
system("wget http://www.regnetworkweb.org/download/human.zip")
system("unzip human.zip")

# Process
TFT_RN<- read.table("human.source",header=FALSE,sep="\t",fill=TRUE,colClasses = "character")
TFT_RN<- TFT_RN[!grepl("hsa-",TFT_RN[,1])&!grepl("hsa-",TFT_RN[,3]),] # Remove miRNA
TFT_RN_ls<-NULL
for(gene in unique(TFT_RN[,1])){
  target_temp<- TFT_RN[TFT_RN[,1]==gene,3]
  TFT_RN_ls<- c(TFT_RN_ls,list(target_temp))
}
names(TFT_RN_ls)<-unique(TFT_RN[,1])
TFT_RN<- TFT_RN_ls

# Some have numbers? 
idx_check<- which(!is.na(as.numeric(names(TFT_RN))))
nr_check<- names(TFT_RN)[idx_check]
nr_check<- nr_check[1:4] # Last one not found 
names_check<- c("GATA1","RARB","VDR","RXRA")
for(i in 1:4){
  name_temp<- names_check[i]
  nr_temp<- nr_check[[i]]
  TFT_RN[[name_temp]]<- union(TFT_RN[[name_temp]],TFT_RN[[nr_temp]])
}
TFT_RN[idx_check]<- NULL # Remove from list

# Save
saveRDS(TFT_RN,file = "TFT_RN.rds")


