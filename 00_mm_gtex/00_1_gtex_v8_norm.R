library(edgeR)

gtex <- read.table("/mnt/lustre/scratch/CBRA/research/projects/heterogeneity_mm/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip=2, head=T)

colnames(gtex) <- gsub("[[:punct:]]","-",colnames(gtex))


### delecte duplicated gene
df_dupl_gen_gtex <- data.frame(table(unlist(lapply(strsplit(paste(gtex[,1]), ".", fixed = T), FUN=function(x) x[1]))))
dupl_gen_gtex <- paste(df_dupl_gen_gtex[df_dupl_gen_gtex$Freq > 1,1])
                                                                                                      
v_menor_var_eli <- c()
for (i in 1:length(dupl_gen_gtex)){
  df_var <- data.frame(apply(gtex[grep(dupl_gen_gtex[i], gtex[,1]),-c(1,2)], 1,var))
  v_menor_var_eli <- c(v_menor_var_eli,rownames(df_var)[which(df_var[,1] ==  min(df_var))])  
}

gtex_hip <- gtex[-as.numeric(v_menor_var_eli),-c(1,2)]
rownames(gtex_hip) <- unlist(lapply(strsplit(gtex[,1][-as.numeric(v_menor_var_eli)], ".", fixed = T), FUN=function(x) x[1])) 


y <- DGEList(counts=gtex_hip)
y <- calcNormFactors(y)
logcpm <- cpm(y, prior.count=3, log=TRUE)

saveRDS(logcpm, file="/mnt/lustre/scratch/CBRA/research/projects/heterogeneity_mm/results/00_mm_gtex/gtex_v8_norm.rds")
