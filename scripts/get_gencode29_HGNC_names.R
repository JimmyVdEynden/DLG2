# Get ENSG genenames fom ht-seq-count output
sampleFiles <- grep("counts",list.files("raw/gene_counts",full.names = T), value=TRUE)
genenames_ENSG<- read.table(sampleFiles[1],stringsAsFactors = F)[,1]

# Add HGNC genenames from biomart
library(biomaRt)
ensembl_hs <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl") # according to https://support.bioconductor.org/p/74906/#74957 
ENSG_HGNC_map_table<- getBM(attributes = c("ensembl_gene_id","ensembl_gene_id_version","external_gene_name","gene_biotype"),filters = "ensembl_gene_id_version",values=genenames_ENSG,mart=ensembl_hs,uniqueRows = TRUE)
rownames(ENSG_HGNC_map_table)<- ENSG_HGNC_map_table[,"ensembl_gene_id_version"]    

# Save
saveRDS(ENSG_HGNC_map_table,file="data/ENSG_HGNC_map_table.rds")
