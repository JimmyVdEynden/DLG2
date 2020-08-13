library(DESeq2)

# Load ENSG_HGNC map
ENSG_HGNC_map_table<- readRDS(file="data/ENSG_HGNC_map_table.rds")

# Load sample information
sample_info<- readRDS(file="data/sample_info.rds")

# Get quantified files
sampleFiles<- grep("counts",list.files("raw/gene_counts",full.names = T), value=TRUE)
sampleNames<- gsub(".*/|\\.gene_counts","",sampleFiles)

# Add condition from sample info and create sample table
sampleCondition <- sample_info[sampleNames,"condition"]
sampleTable<-data.frame(sampleName=sampleNames, fileName=sampleFiles, condition=sampleCondition)

# From this fase combina with sample information!
ddsHTSeq<- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, design=~condition)
colData(ddsHTSeq)$condition<- factor(colData(ddsHTSeq)$condition, levels=c("Ctrl","Ctrl_RA","DLG2","DLG2_RA"))
dds<- DESeq(ddsHTSeq)

# Diff expression analysis for each contion
res_RA_Ctrl<- results(dds, contrast = c("condition","Ctrl_RA","Ctrl"))
res_DLG2_Ctrl<- results(dds, contrast = c("condition","DLG2","Ctrl"))
res_DLG2RA_Ctrl<- results(dds, contrast = c("condition","DLG2_RA","Ctrl"))
res_DLG2RA_CtrlRA<- results(dds, contrast = c("condition","DLG2_RA","Ctrl_RA"))
res_DLG2RA_DLG2<- results(dds, contrast = c("condition","DLG2_RA","DLG2"))

# Process each conditions 1) 
results_diff_expr<- list()
for(i in 1:5){
  res<- list(res_RA_Ctrl,res_DLG2_Ctrl,res_DLG2RA_Ctrl,res_DLG2RA_CtrlRA,res_DLG2RA_DLG2)[[i]]
  comparison<- c("Ctrl+RA Vs. Ctrl","DLG2 Vs. Ctrl","DLG2+RA Vs. Ctrl","DLG2+RA Vs. Ctrl+RA","DLG2+RA Vs. DLG2")[i]
  
  # Only coding genes
  res_proc<- res[ENSG_HGNC_map_table[rownames(res),"gene_biotype"]=="protein_coding"&!is.na(ENSG_HGNC_map_table[rownames(res),"gene_biotype"]),]
  res_proc<- res[!is.na(ENSG_HGNC_map_table[rownames(res),"external_gene_name"]),]
  
  # Convert to HGNC names
  # No one on one mapping!!!! Multiple gene IDs --> Take the one that has (highest) expression data
  res_proc_HGNC_names<- ENSG_HGNC_map_table[rownames(res_proc),"external_gene_name"]
  res_proc$HGNC<- res_proc_HGNC_names
  res_proc<- res_proc[order(res_proc$HGNC,res_proc[,"baseMean"],decreasing=T),]
  res_proc<- res_proc[!duplicated(res_proc$HGNC),]
  res_proc$ENSG<- rownames(res_proc)
  rownames(res_proc)<- res_proc$HGNC
  
  # Sort on genenames
  res_proc<- res_proc[order(res_proc[,"HGNC"]),]

  # Redo padj for selected (coding) genes
  # res_proc$padj<- p.adjust(res_proc$pvalue,"fdr")
  
  # Save to dataframe
  assign(paste0("results",i),as.data.frame(res_proc)) 
  results_diff_expr<- c(results_diff_expr,res_proc)
}

# Get normalized counts for each sample & save
normalized_counts <- counts(dds, normalized=T)
normalized_counts<- normalized_counts[res_proc$ENSG,]
rownames(normalized_counts)<- res_proc$HGNC
normalized_counts<- as.data.frame(normalized_counts)
WriteXLS::WriteXLS("normalized_counts","results/tables/DLG2_normalized_counts.xlsx",row.names = T)
saveRDS(normalized_counts,file = "data/DLG2_normalized_counts.rds")

# Save results
results_all<- c("results1","results2","results3","results4","results5")
results_all_sheetnames<-  c("Ctrl+RA Vs. Ctrl","DLG2 Vs. Ctrl","DLG2+RA Vs. Ctrl","DLG2+RA Vs. Ctrl+RA","DLG2+RA Vs. DLG2") 
WriteXLS::WriteXLS(results_all, "results/tables/DLG2_diff_expression_results.xlsx",SheetNames = results_all_sheetnames,row.names = T)

names(results_diff_expr)<- results_all_sheetnames
saveRDS(results_diff_expr,file = "data/DLG2_diff_expression_results.rds")
                   


