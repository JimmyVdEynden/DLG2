#################################
# proteomics_process.R
#################################

# Save quantified data as df
#############################

DLG2_prot<- readxl::read_excel("raw/proteomics/3062_Siaw_QMS_Fusion_190123_06-25-(1)_1%FDR_proteins_190206.xlsx",skip = 2)
DLG2_prot_df<- as.data.frame(DLG2_prot)

# Extract HGNC IDs and add as rownames
HGNC_ids<- DLG2_prot_df[,"Description"]
HGNC_ids<- gsub(".*GN=","",HGNC_ids)
HGNC_ids<- gsub(" PE=.*","",HGNC_ids)
DLG2_prot_df<- cbind(DLG2_prot_df,HGNC_ids)
DLG2_prot_df<- DLG2_prot_df[!duplicated(DLG2_prot_df[,"HGNC_ids"]),] # Only 10 duplicated
rownames(DLG2_prot_df)<- DLG2_prot_df[,"HGNC_ids"]

# Only extract normalized data  
DLG2_prot_df<- DLG2_prot_df[,grepl("(Normalized)",colnames(DLG2_prot))]
colnames(DLG2_prot_df)<- gsub(".*, ","",colnames(DLG2_prot_df))
colnames(DLG2_prot_df)<- gsub("GFP-DMSO ","Ctrl_",colnames(DLG2_prot_df))
colnames(DLG2_prot_df)<- gsub("GFP-RA ","RA_",colnames(DLG2_prot_df))
colnames(DLG2_prot_df)<- gsub("DLG2-DMSO ","DLG2_",colnames(DLG2_prot_df))
colnames(DLG2_prot_df)<- gsub("DLG2-RA ","DLG2_RA_",colnames(DLG2_prot_df))

#Save
saveRDS(DLG2_prot_df,file = "data/DLG2_prot.rds")

# DE
####
library(ROTS)

# Create for all conditions
results_diff_expr<- list()
for(i in 1:5){
  cat(i, " ")
  cond<- c("RA","DLG2","DLG2_RA","DLG2_RA","DLG2_RA")[i]
  ctrl<- c("Ctrl","Ctrl","Ctrl","RA","DLG2")[i]
  if(cond=="DLG2_RA") n_cond<-2
  else n_cond<- 3
  input<- DLG2_prot_df[,c(paste(cond,1:n_cond,sep="_"),paste(ctrl,1:3,sep="_"))]
  groups<- as.numeric(factor(gsub("_[[:digit:]]","",colnames(input)),levels = c(cond,ctrl)))
  input<- input[!rownames(input)%in%c("KRT83","KRT33A"),] # too many missing values
  results<- ROTS(data = input, groups = groups , B = 1000 , K = 500 , seed = 1234,log=F) # Put log to false!!!
  results_df<- summary(results,fdr = 1)
  results_df<- cbind(results_df,logFC=results$logfc[results_df[,"Row"]])
  assign(paste("results",ctrl,cond,sep="_"),as.data.frame(results_df)) 
  results_diff_expr<- c(results_diff_expr,list(results_df))
}

# Save results
cond<- c("RA","DLG2","DLG2_RA","DLG2_RA","DLG2_RA")
ctrl<- c("Ctrl","Ctrl","Ctrl","RA","DLG2")
results_all<- paste("results",ctrl,cond,sep="_")
results_all_sheetnames<-  c("Ctrl+RA Vs. Ctrl","DLG2 Vs. Ctrl","DLG2+RA Vs. Ctrl","DLG2+RA Vs. Ctrl+RA","DLG2+RA Vs. DLG2") 
WriteXLS::WriteXLS(results_all, "results/tables/DLG2_diff_expression_proteomics_results.xlsx",SheetNames = results_all_sheetnames,row.names = T)

names(results_diff_expr)<- results_all_sheetnames
saveRDS(results_diff_expr,file = "data/DLG2_diff_expression_proteomics_results.rds")
