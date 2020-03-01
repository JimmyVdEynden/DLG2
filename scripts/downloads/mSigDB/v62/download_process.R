####################
# download_process.R
####################

# Download latest manually
# GO database from MSigDB at http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2/c5.all.v6.2.symbols.gmt
# Reactome database: http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2/c2.cp.reactome.v6.2.symbols.gmt
# Kegg: http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2/c2.cp.kegg.v6.2.symbols.gmt

# Read data
GO<-read.table("c5.all.v6.2.symbols.gmt",header=FALSE,sep="\t",colClasses = rep("character",2133),row.names = 1,fill=TRUE)
Rea<-read.table("c2.cp.reactome.v6.2.symbols.gmt",header=FALSE,sep="\t",colClasses = rep("character",2133),row.names = 1,fill=TRUE)
Kegg<-read.table("c2.cp.kegg.v6.2.symbols.gmt",header=FALSE,sep="\t",colClasses = rep("character",2133),row.names = 1,fill=TRUE)

# Make lists
GO_ls<- list()
for(i in 1:nrow(GO)){
  cat(i," ")
  GO_name<- rownames(GO)[i]
  GO_ls[[GO_name]]<- setdiff(unique(as.character(GO[GO_name,-1])),"") 
}

Rea_ls<- list()
for(i in 1:nrow(Rea)){
  cat(i," ")
  Rea_name<- rownames(Rea)[i]
  Rea_ls[[Rea_name]]<- setdiff(unique(as.character(Rea[Rea_name,-1])),"") 
}

Kegg_ls<- list()
for(i in 1:nrow(Kegg)){
  cat(i," ")
  Kegg_name<- rownames(Kegg)[i]
  Kegg_ls[[Kegg_name]]<- setdiff(unique(as.character(Kegg[Kegg_name,-1])),"") 
}

# Save
save(GO_ls,Kegg_ls,Rea_ls,file = "MSigDB_v62.RData")

# Add 3 types of GO:
# bp= http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2/c5.bp.v6.2.symbols.gmt
# cc= http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2/c5.cc.v6.2.symbols.gmt
# mf= http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2/c5.mf.v6.2.symbols.gmt
GO_cc<-read.table("c5.cc.v6.2.symbols.gmt",header=FALSE,sep="\t",colClasses = rep("character",2133),row.names = 1,fill=TRUE)
GO_bp<-read.table("c5.bp.v6.2.symbols.gmt",header=FALSE,sep="\t",colClasses = rep("character",2133),row.names = 1,fill=TRUE)
GO_mf<-read.table("c5.mf.v6.2.symbols.gmt",header=FALSE,sep="\t",colClasses = rep("character",2133),row.names = 1,fill=TRUE)

# Make lists
GO_cc_ls<- list()
for(i in 1:nrow(GO_cc)){
  cat(i," ")
  GO_name<- rownames(GO_cc)[i]
  GO_cc_ls[[GO_name]]<- setdiff(unique(as.character(GO_cc[GO_name,-1])),"") 
}

GO_bp_ls<- list()
for(i in 1:nrow(GO_bp)){
  cat(i," ")
  GO_name<- rownames(GO_bp)[i]
  GO_bp_ls[[GO_name]]<- setdiff(unique(as.character(GO_bp[GO_name,-1])),"") 
}

GO_mf_ls<- list()
for(i in 1:nrow(GO_mf)){
  cat(i," ")
  GO_name<- rownames(GO_mf)[i]
  GO_mf_ls[[GO_name]]<- setdiff(unique(as.character(GO_mf[GO_name,-1])),"") 
}

save(GO_cc_ls,GO_bp_ls,GO_mf_ls,file = "MSigDB_GO_3types_v62.RData")

