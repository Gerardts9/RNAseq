library(data.table);library(openxlsx);library(edgeR);library(preprocessCore);library(qvalue);library(tidyverse)


### Plot residuals:

ann <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/Annotation/gencode.v26.annotation.fixed.gtf.gz")

##### By Case:

########################################
counts.tpm <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/Gene_counts/Counts.RSEM.All.txt")

tec <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/Tables/Technical_Complete.xlsx")
tec <- tec[tec$Seq_ID %in% colnames(counts.tpm),]

bio <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/20221115_Datos Pacientes.xlsx")
bio <- bio[2:nrow(bio),]
bio <- bio[bio$Muestra %in% colnames(counts.tpm),]

table(bio$Status)

# Smoking:
bio$Smoking2 <- ifelse(bio$Smoking == 0, 0, 1)


dge <- DGEList(counts = counts.tpm)

# Set "Control" to be the first level:
AAA <- factor(x = c(rep("Case", 96), rep("Control", 44)),
              levels = c("Control", "Case"))

keep <- apply(dge$counts, 1, function(x) sum(x > 0.5) >= 70)
dge <- dge[keep,,keep.lib.sizes=FALSE]

paste0("We keep ", sum(keep), " genes") # 14,675
paste0("We remove ", sum(!keep), " genes") # 12,615

counts.tpm <- as.data.frame(normalize.quantiles(as.matrix(dge$counts), keep.names = T))

tpm_matrix <- as.matrix(counts.tpm)

# Retrieve p-values for Case/Control variable:
Pvalue <-
  apply(counts.tpm, 1, function(x)
    summary(
      lm(
        x ~ tec$Case + tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO
      )
    )[4][[1]][2, 4])


# Retrieve betas for Case/Control variable:
Coef <-
  apply(counts.tpm, 1, function(x)
    summary(
      lm(
        x ~ tec$Case + tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO
      )
    )[4][[1]][2, 1])

# Create q-value object:
qobj <- qvalue(p = Pvalue)

Results <- data.frame(
  pvalue = Pvalue,
  qvalue = qobj$qvalues,
  fdr = p.adjust(Pvalue, method = "fdr"),
  pi1 = 1-qobj$pi0,
  Beta = Coef
)

Results <- merge(Results, ann[,c("gene_id","gene_name")], by.x = 0, by.y = "gene_id")
names(Results)[names(Results) == "Row.names"] <- "gene_id"
Results <- Results[!duplicated(Results$gene_name),]
Results <- Results %>% remove_rownames %>% column_to_rownames(var="gene_name")
Results <- Results[order(Results$pvalue),]

nrow(Results[Results$fdr < 0.05,]) # 7454

# Use gene name instead of ENSEMBL:

counts <- merge(counts.tpm, ann[,c("gene_id","gene_name")], by.x = 0, by.y = "gene_id")
counts <- counts[,-1]
counts <- counts[!duplicated(counts$gene_name),]
counts <- counts %>% remove_rownames %>% column_to_rownames(var="gene_name")

counts <- as.matrix(counts)

cases.tpm <- counts[,1:96]
controls.tpm <- counts[,97:ncol(counts)]

dim(cases.tpm)
dim(controls.tpm)

table(AAA)

glist <- c("EXTL3","ZFR","VPS37C","ZNF784","DUSP8","USP33","RFX1","DISP1")


for (gene in glist) {
  try(res <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res <-  data.frame(Expression = res)
  
  y_coord <- max(res$Expression)
  
  text <- format(as.numeric(Results[rownames(Results) %in% gene,]$fdr), digits = 2)
  p <- ggplot(res, aes(x = AAA, y = Expression, fill = AAA)) + 
    geom_boxplot(outlier.shape = NA) + 
    ggtitle(gene) +
    ylab("Expression (Normalized TPM)") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 22),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20)
    ) +
    scale_fill_manual(values=c("cyan3", "red")) +
    annotate("text", x = 1.5, y = y_coord, label = paste0("FDR P-value:", text)) +
    geom_jitter(data = res, aes(x = AAA, y = Expression), color = "black", size = 1, alpha = 0.5)
  
  p
  
  ggsave(paste0("C://Users/Gerard/Desktop/", gene, ".png"), width = 7, height = 5)
}
########################################



#### By diameter:
########################################
counts.tpm <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/Gene_counts/Counts.RSEM.Cases.txt")

tec <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/Tables/First Analysis/Technical_CompleteV2.xlsx")
tec <- tec[tec$SampleID %in% colnames(counts.tpm),]

bio <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/20221115_Datos Pacientes.xlsx")
bio <- bio[bio$Muestra %in% colnames(counts.tpm),]

# Smoking:
bio$Smoking2 <- ifelse(bio$Smoking == 0, 0, 1)

dge <- DGEList(counts = counts.tpm)

keep <- apply(dge$counts, 1, function(x) sum(x > 0.5) >= 48)
dge <- dge[keep,,keep.lib.sizes=FALSE]

paste0("We keep ", sum(keep), " genes") # 14,779
paste0("We remove ", sum(!keep), " genes") # 12,511


counts.tpm <- as.data.frame(normalize.quantiles(as.matrix(dge$counts), keep.names = TRUE))


# Retrieve p-values for Diameter variable:
Pvalue <-
  apply(counts.tpm, 1, function(x)
    summary(
      lm(
        x ~ bio$aortc_diameter_mm + tec$Date + tec$Batch + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO
      )
    )[4][[1]][2, 4])

Coef <-
  apply(counts.tpm, 1, function(x)
    summary(
      lm(
        x ~ bio$aortc_diameter_mm + tec$Date + tec$Batch + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO
      )
    )[4][[1]][2, 1])

qobj <- qvalue(p = Pvalue)

Results <- data.frame(
  pvalue = Pvalue,
  qvalue = qobj$qvalues,
  fdr = p.adjust(Pvalue, method = "fdr"),
  pi1 = 1-qobj$pi0,
  Beta = Coef
)

Results <- merge(Results, ann[,c("gene_id","gene_name")], by.x = 0, by.y = "gene_id")
names(Results)[names(Results) == "Row.names"] <- "gene_id"
Results <- Results[!duplicated(Results$gene_name),]
Results <- Results %>% remove_rownames %>% column_to_rownames(var="gene_name")
Results <- Results[order(Results$pvalue),]

nrow(Results[Results$fdr < 0.05,]) # 32


counts <- merge(counts.tpm, ann[,c("gene_id","gene_name")], by.x = 0, by.y = "gene_id")
counts <- counts[,-1]
counts <- counts[!duplicated(counts$gene_name),]
counts <- counts %>% remove_rownames %>% column_to_rownames(var="gene_name")

counts <- as.matrix(counts)

dia <- bio[!is.na(bio$aortc_diameter_mm),][,c("Muestra","aortc_diameter_mm")]

counts <- counts[,colnames(counts) %in% dia$Muestra]
tec2 <- tec[tec$SampleID %in% colnames(counts),]
bio2 <- bio[bio$Muestra %in% colnames(counts),]

head(counts)

dia <- dia %>% remove_rownames %>% column_to_rownames(var="Muestra")

glist <- c("EXTL3","ZFR","VPS37C","ZNF784","DUSP8","USP33","RFX1","DISP1")

for (gene in glist) {
  res <- resid(lm(counts[gene,] ~ tec2$Date + tec2$Batch + tec2$GC_Mean + tec2$RIN + tec2$DV200 + tec2$Qubit + bio2$age + bio2$SEXO))
  res <- data.frame(Expression = res,
                    Diameter = dia)

  text <- format(as.numeric(Results[rownames(Results) %in% gene,]$fdr), digits = 2)

  y_coord <- max(res$Expression)
  x_coord <- mean(res$aortc_diameter_mm)

  ggplot(res, aes(x = aortc_diameter_mm, y = Expression)) + geom_point(color = "orangered2") + ggtitle(paste(gene)) + ylab ("Expression (Normalized TPM)") + xlab("Aortic Diameter (mm)") +
    geom_smooth(method = 'lm', formula = y ~ x, color = "black") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 22),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")) +
    annotate("text", x = x_coord, y = y_coord, label = paste0("FDR P-value:", text))

  ggsave(paste0("C://Users/Gerard/Desktop/AAA/RNAseq/Plots/", gene, "_Diameter.png"), width = 7, height = 5)
}
########################################




### Daniel Genes:


#glist <- c("PDP1","PDP2","PDHA1","PDK1","PDK2","PDK3","PDK4","ACO1","ACO2","IDH1","IDH2","IDH3A","IDH3B","OGDH","SUCLA2",
#"SDHA","SDHB","SDHC","SDHD","GLS","GLUL","GLUD1","GAD1","GAD2","ABAT","ALDH5A1")



