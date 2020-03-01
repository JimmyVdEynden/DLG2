# Get information manually from 2 excels sent by Joachim and Ganesh on 2/11/18:
# see raw/sample_info/

sample_name<- c('GFP_5Y1',
                'GFP_5Y2',
                'GFP_5Y3',
                'GFP_5Y4',
                'GFP_5Y5',
                'GFP_5Y6',
                'GFP_5Y7',
                'GFP_5Y8',
                'DLG2_5Y1',
                'DLG2_5Y2',
                'DLG2_5Y3',
                'DLG2_5Y4',
                'DLG2_5Y5',
                'DLG2_5Y6',
                'DLG2_5Y7',
                'DLG2_5Y8'
                )
hasDLG2Overexpr<- c(rep(F,8),rep(T,8)) 
isRA<- c(rep(F,4),rep(T,4),rep(F,4),rep(T,4)) 
condition<- c(rep("Ctrl",4),rep("Ctrl_RA",4),rep("DLG2",4),rep("DLG2_RA",4))
sample_matrix<- data.frame(hasDLG2Overexpr,isRA,condition)
rownames(sample_matrix)<- sample_name
saveRDS(sample_matrix,file="data/sample_info.rds")
