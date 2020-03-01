# Download markers from nature genetics paper 
# wget https://media.nature.com/original/nature-assets/ng/journal/v49/n8/extref/ng.3899-S3.xlsx

# Process
# library(gdata)
# diff_marker<- read.xls("downloads/diff_markers/ng.3899-S3.xlsx")
diff_marker<- as.data.frame(readxl::read_xlsx("ng.3899-S3.xlsx"))
MES<- as.character(diff_marker[!is.na(diff_marker[,2])&diff_marker[,2]=="MES",1])
ADRN<- as.character(diff_marker[!is.na(diff_marker[,2])&diff_marker[,2]=="ADRN",1])

# Save
save(MES,ADRN,file = "vGron_diff_markers.RData")
