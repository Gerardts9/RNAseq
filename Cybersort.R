library(data.table);library(edgeR);library(tidyverse);library(openxlsx);library(preprocessCore);library("org.Hs.eg.db");library(stats);library(qvalue);library(ggpubr);library(VennDiagram)
library(readxl)



lm22 <- data.frame(read_xls("C://Users/Gerard/Downloads/41592_2015_BFnmeth3337_MOESM207_ESM.xls", sheet = 1, skip = 13))

head(lm22)

fwrite(lm22[,3:ncol(lm22)], "C://Users/Gerard/Desktop/AAA/RNAseq/Cibersort/lm22.txt", sep = "\t")




ann <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/Annotation/gencode.v26.annotation.fixed.gtf.gz")
counts.tpm <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/Gene_counts/Counts.RSEM.All.txt")
tpm_matrix <- as.matrix(counts.tpm)


counts.tpm <- merge(counts.tpm, ann[,c("gene_id","gene_name")], by.x = 0, by.y = "gene_id")
counts.tpm <- counts.tpm[,-1]
counts.tpm <- counts.tpm[!duplicated(counts.tpm$gene_name),]
counts.tpm <- counts.tpm %>% remove_rownames %>% column_to_rownames(var="gene_name")

head(counts.tpm)

fwrite(counts.tpm, "C://Users/Gerard/Desktop/AAA/RNAseq/Cibersort/counts.txt", sep = "\t", row.names = T, col.names = T)



tec <- read.xlsx("/home/gerard/AAA/Tables/Technical_Complete.xlsx")
tec <- tec[tec$Seq_ID %in% colnames(counts.tpm),]

bio <- read.xlsx("/home/gerard/AAA/Tables/20221115_Datos Pacientes.xlsx")
bio <- bio[2:nrow(bio),]
bio <- bio[bio$Muestra %in% colnames(counts.tpm),]

table(bio$Status)

# Smoking:
bio$Smoking2 <- ifelse(bio$Smoking == 0, 0, 1)

dge <- DGEList(counts = counts.tpm)

# Set "Control" to be the first level:
AAA <- factor(x = c(rep("Case", 96), rep("Control", 44)),
              levels = c("Control", "Case"))

# Batch and FC are not included because are highly correlated with AAA:
design <-
  model.matrix(~ AAA + tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$SEXO + bio$age + tec$Library)

colnames(design) <-
  c("(Intercept)",
    "AAA",
    "FC1",
    "FC2",
    "Lane1",
    "Lane2",
    "GC_Mean",
    "RIN",
    "DV200",
    "Qubit",
    "Sex",
    "Age",
    "Library")

keep <- apply(dge$counts, 1, function(x) sum(x > 0.5) >= 70)
dge <- dge[keep,,keep.lib.sizes=FALSE]

paste0("We keep ", sum(keep), " genes") # 14,675
paste0("We remove ", sum(!keep), " genes") # 12,615 

counts.tpm <- as.data.frame(normalize.quantiles(as.matrix(dge$counts), keep.names = T))

tpm_matrix <- as.matrix(counts.tpm)
