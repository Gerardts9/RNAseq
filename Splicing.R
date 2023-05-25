library(data.table);library(clusterProfiler);library(openxlsx)

med <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 1)
adv <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 2)

ann <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/Annotation/gencode.v26.annotation.fixed.gtf.gz")
ann$gene_id2 <- substr(ann$gene_id, 1, 15)

res <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/SUPPA/AAA_Results.dpsi")

head(res)

nrow(res) # 179,651

nrow(na.omit(res)) # 45,577

res2 <- na.omit(res)

names(res2)

res2$fdr <- p.adjust(res2$AAA_Cases.AAA_Controls_p.val, method = "fdr")
nrow(res2[res2$fdr < 0.05,])

sig <- res2[abs(res2$AAA_Cases.AAA_Controls_dPSI) > 0.1 & res2$AAA_Cases.AAA_Controls_p.val < 0.05, ]
nrow(sig)


#write.table(sig, "/home/gerard/AAA/suppa/Results/AAA_Sig_Results.txt") # 173

sig$event <- rownames(sig)
sig$gene_id2 <- substr(rownames(sig), 1, 15)
sig <- merge(sig, ann, by = "gene_id2")

sig$gene_name

cc <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Cases_Controls.txt")
cc$gene_id2 <- substr(cc$gene_id, 1, 15)

length(intersect(sig$gene_id2, cc$gene_id2))/length(unique(sig$gene_id2)) # 50% Of the significant splicing are DEG.

intersect(sig$gene_name, cc$V1)

sig <- sig[order(sig$AAA_Cases.AAA_Controls_p.val),]



