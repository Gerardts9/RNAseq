library(data.table);library(openxlsx);library(edgeR);library(preprocessCore);library(qvalue);library(tidyverse);library(ggsignif)
library(RColorBrewer);library(matrixStats);library(pheatmap);library(pheatmap);library(RColorBrewer)


### Plot residuals:

ann <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/Annotation/gencode.v26.annotation.fixed.gtf.gz")
head(ann)

# Define function to detect if number is in scientific notation:
is_scientific <- function(x) {
  return(grepl("[eE]", x))
}


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

dge2 <- DGEList(counts = counts.tpm)
dge.rm <- dge2[!keep,,keep.lib.sizes=FALSE]

paste0("We keep ", sum(keep), " genes") # 14,675
paste0("We remove ", sum(!keep), " genes") # 12,615

counts.tpm <- as.data.frame(normalize.quantiles(as.matrix(dge$counts), keep.names = T))
counts.tpm.rm <- as.data.frame(normalize.quantiles(as.matrix(dge.rm$counts), keep.names = T))

dim(counts.tpm)
dim(counts.tpm.rm)


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

nrow(Results[Results$fdr < 0.05,]) # 7,454

# Use gene name instead of ENSEMBL for lowely expressed genes:
counts.rm <- merge(counts.tpm.rm, ann[,c("gene_id","gene_name")], by.x = 0, by.y = "gene_id")
counts.rm <- counts.rm[,-1]
counts.rm <- counts.rm[!duplicated(counts.rm$gene_name),]
counts.rm <- counts.rm %>% remove_rownames %>% column_to_rownames(var="gene_name")


# Use gene name instead of ENSEMBL for significant genes:
counts <- merge(counts.tpm, ann[,c("gene_id","gene_name")], by.x = 0, by.y = "gene_id")
counts <- counts[,-1]
counts <- counts[!duplicated(counts$gene_name),]
counts <- counts %>% remove_rownames %>% column_to_rownames(var="gene_name")

counts <- as.matrix(counts)
dim(counts)
dim(counts.tpm)

cases.tpm <- counts[,1:96]
controls.tpm <- counts[,97:ncol(counts)]

dim(cases.tpm)
dim(controls.tpm)

table(AAA)
########################################



# 1) Comparison plots between AAA and control individuals:

# 1.1) Mitophagy Related Genes (Park2 changed to PARK2) (Prohibitn2 changed to PHB2):
########################################
glist <- c(
  "PARK2", "PINK1", "USP30", "BNIP3", "BNIP3L", "FUNDC1", "BECN1", "AMBRA1",
  "OPTN", "PHB2", "PPARGC1A", "ATG5", "ATG7", "MAP1LC3B"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}


glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)

for (gene in glist) {
  p_value <- Results[gene, ]$fdr
  
  if (is_scientific(p_value)) {
    annots <- signif(p_value, digits = 3)
  }
  else {
    annots <- round(p_value, 4)
  }
  
  ggplot(res, aes(x = AAA, y = res[, gene], fill = AAA)) +
    geom_boxplot() +
    ggtitle(gene) +
    ylab("Expression (Normalized TPM)") +
    labs(x = gene) +
    theme(
      plot.title = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.margin = margin(t = 0),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.title.x = element_text(
        size = 12,
        angle = 45,
        vjust = 0.5
      ),
      panel.background = element_blank()
    ) +
    scale_fill_manual(values = c("turquoise2", "salmon2"),
                      labels = c("Controls", "AAA")) +
    geom_signif(
      comparisons = list(c("Control", "Case")),
      map_signif_level = TRUE,
      textsize = 4,
      test = "t.test",
      annotations = annots
    )
  ggsave(
    paste0(
      "C://Users/Gerard/Desktop/Plots Hui/Results/Plots_Mitophagy/",
      gene,
      "_Plot.png"
    ),
    width = 7,
    height = 5
  )
}
########################################



# 1.2) Mitochondrial Function Related Genes:

# 1.2.1) Complex I (Presence of MT genes, MT genes are level 3, and removed consequently):
########################################
glist <- c(
  "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7", "NDUFA8", 
  "NDUFA9", "NDUFA10", "NDUFA11", "NDUFA12", "NDUFA13", "NDUFB1", "NDUFB2", 
  "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", 
  "NDUFB11", "NDUFC1", "NDUFC2", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", 
  "NDUFS6", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2", "NDUFV3", "MT-ND1", "MT-ND2", 
  "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}

glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)
table(res$Group)

for (gene in glist) {
  p_value <- Results[gene, ]$fdr
  
  if (is_scientific(p_value)) {
    annots <- signif(p_value, digits = 3)
  }
  else {
    annots <- round(p_value, 4)
  }
  
  ggplot(res, aes(x = AAA, y = res[, gene], fill = AAA)) +
    geom_boxplot() +
    ggtitle(gene) +
    ylab("Expression (Normalized TPM)") +
    labs(x = gene) +
    theme(
      plot.title = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.margin = margin(t = 0),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.title.x = element_text(
        size = 12,
        angle = 45,
        vjust = 0.5
      ),
      panel.background = element_blank()
    ) +
    scale_fill_manual(values = c("turquoise2", "salmon2"),
                      labels = c("Controls", "AAA")) +
    geom_signif(
      comparisons = list(c("Control", "Case")),
      map_signif_level = TRUE,
      textsize = 4,
      test = "t.test",
      annotations = annots
    )
  ggsave(
    paste0(
      "C://Users/Gerard/Desktop/Plots Hui/Results/Plots_Mithocondrial_Complex_I/",
      gene,
      "_Plot.png"
    ),
    width = 7,
    height = 5
  )
}
########################################


# 1.2.2) Complex II:
########################################
glist <- c(
  "SDHA","SDHB","SDHC","SDHD"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}

glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)
table(res$Group)

for (gene in glist) {
  p_value <- Results[gene, ]$fdr
  
  if (is_scientific(p_value)) {
    annots <- signif(p_value, digits = 3)
  }
  else {
    annots <- round(p_value, 4)
  }
  
  ggplot(res, aes(x = AAA, y = res[, gene], fill = AAA)) +
    geom_boxplot() +
    ggtitle(gene) +
    ylab("Expression (Normalized TPM)") +
    labs(x = gene) +
    theme(
      plot.title = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.margin = margin(t = 0),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.title.x = element_text(
        size = 12,
        angle = 45,
        vjust = 0.5
      ),
      panel.background = element_blank()
    ) +
    scale_fill_manual(values = c("turquoise2", "salmon2"),
                      labels = c("Controls", "AAA")) +
    geom_signif(
      comparisons = list(c("Control", "Case")),
      map_signif_level = TRUE,
      textsize = 4,
      test = "t.test",
      annotations = annots
    )
  ggsave(
    paste0(
      "C://Users/Gerard/Desktop/Plots Hui/Results/Plots_Mithocondrial_Complex_II/",
      gene,
      "_Plot.png"
    ),
    width = 7,
    height = 5
  )
}
########################################



# 1.2.3) Complex III (Presence of ONE MT gene (level 3)):
########################################
glist <- c(
  "UQCRC1", "UQCRC2", "UQCRB", "UQCRQ", "UQCRH", "UQCRFS1", "MT-CYB"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}

glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)
table(res$Group)

for (gene in glist) {
  p_value <- Results[gene, ]$fdr
  
  if (is_scientific(p_value)) {
    annots <- signif(p_value, digits = 3)
  }
  else {
    annots <- round(p_value, 4)
  }
  
  ggplot(res, aes(x = AAA, y = res[, gene], fill = AAA)) +
    geom_boxplot() +
    ggtitle(gene) +
    ylab("Expression (Normalized TPM)") +
    labs(x = gene) +
    theme(
      plot.title = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.margin = margin(t = 0),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.title.x = element_text(
        size = 12,
        angle = 45,
        vjust = 0.5
      ),
      panel.background = element_blank()
    ) +
    scale_fill_manual(values = c("turquoise2", "salmon2"),
                      labels = c("Controls", "AAA")) +
    geom_signif(
      comparisons = list(c("Control", "Case")),
      map_signif_level = TRUE,
      textsize = 4,
      test = "t.test",
      annotations = annots
    )
  ggsave(
    paste0(
      "C://Users/Gerard/Desktop/Plots Hui/Results/Plots_Mithocondrial_Complex_III/",
      gene,
      "_Plot.png"
    ),
    width = 7,
    height = 5
  )
}
########################################



# 1.2.4) Complex IV (Presence of MT genes (level 3)) (COX6A2 and COX6B2 removed for low expression):
########################################
glist <- c(
  "COX4I1", "COX4I2", "COX5A", "COX5B", "COX6A1", "COX6A2", "COX6B1", "COX6B2", 
  "COX6C", "COX7A1", "COX7A2", "COX7A2L", "COX7B", "COX7C", "COX8A", "MT-CO1", 
  "MT-CO2", "MT-CO3"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}

glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)
table(res$Group)

for (gene in glist) {
  p_value <- Results[gene, ]$fdr
  
  if (is_scientific(p_value)) {
    annots <- signif(p_value, digits = 3)
  }
  else {
    annots <- round(p_value, 4)
  }
  
  ggplot(res, aes(x = AAA, y = res[, gene], fill = AAA)) +
    geom_boxplot() +
    ggtitle(gene) +
    ylab("Expression (Normalized TPM)") +
    labs(x = gene) +
    theme(
      plot.title = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.margin = margin(t = 0),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.title.x = element_text(
        size = 12,
        angle = 45,
        vjust = 0.5
      ),
      panel.background = element_blank()
    ) +
    scale_fill_manual(values = c("turquoise2", "salmon2"),
                      labels = c("Controls", "AAA")) +
    geom_signif(
      comparisons = list(c("Control", "Case")),
      map_signif_level = TRUE,
      textsize = 4,
      test = "t.test",
      annotations = annots
    )
  ggsave(
    paste0(
      "C://Users/Gerard/Desktop/Plots Hui/Results/Plots_Mithocondrial_Complex_IV/",
      gene,
      "_Plot.png"
    ),
    width = 7,
    height = 5
  )
}
########################################



# 1.2.5) Complex V (changed all names to the annotation ones):
########################################
glist <- c(
  "ATP5A1", "ATP5B", "ATP5C1", "ATP5D", "ATP5E", "ATP5G1", "ATP5G2", 
  "ATP5G3", "ATP5I", "ATP5J2", "ATP5L", "ATP5F1", "ATP5H", "ATP5J", 
  "MT-ATP6", "MT-ATP8"
)


# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}

glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)
table(res$Group)

for (gene in glist) {
  p_value <- Results[gene, ]$fdr
  
  if (is_scientific(p_value)) {
    annots <- signif(p_value, digits = 3)
  }
  else {
    annots <- round(p_value, 4)
  }
  
  ggplot(res, aes(x = AAA, y = res[, gene], fill = AAA)) +
    geom_boxplot() +
    ggtitle(gene) +
    ylab("Expression (Normalized TPM)") +
    labs(x = gene) +
    theme(
      plot.title = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.margin = margin(t = 0),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.title.x = element_text(
        size = 12,
        angle = 45,
        vjust = 0.5
      ),
      panel.background = element_blank()
    ) +
    scale_fill_manual(values = c("turquoise2", "salmon2"),
                      labels = c("Controls", "AAA")) +
    geom_signif(
      comparisons = list(c("Control", "Case")),
      map_signif_level = TRUE,
      textsize = 4,
      test = "t.test",
      annotations = annots
    )
  ggsave(
    paste0(
      "C://Users/Gerard/Desktop/Plots Hui/Results/Plots_Mithocondrial_Complex_V/",
      gene,
      "_Plot.png"
    ),
    width = 7,
    height = 5
  )
}
########################################


# 1.2.6) Others (ATP5PF removed because repeated in previous set of genes):
########################################
glist <- c(
  "CYCS", "NDUFAF1", "SCO1", "SCO2"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}

glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)
table(res$Group)

for (gene in glist) {
  p_value <- Results[gene, ]$fdr
  
  if (is_scientific(p_value)) {
    annots <- signif(p_value, digits = 3)
  }
  else {
    annots <- round(p_value, 4)
  }
  
  ggplot(res, aes(x = AAA, y = res[, gene], fill = AAA)) +
    geom_boxplot() +
    ggtitle(gene) +
    ylab("Expression (Normalized TPM)") +
    labs(x = gene) +
    theme(
      plot.title = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.margin = margin(t = 0),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.title.x = element_text(
        size = 12,
        angle = 45,
        vjust = 0.5
      ),
      panel.background = element_blank()
    ) +
    scale_fill_manual(values = c("turquoise2", "salmon2"),
                      labels = c("Controls", "AAA")) +
    geom_signif(
      comparisons = list(c("Control", "Case")),
      map_signif_level = TRUE,
      textsize = 4,
      test = "t.test",
      annotations = annots
    )
  ggsave(
    paste0(
      "C://Users/Gerard/Desktop/Plots Hui/Results/Plots_Mithocondrial_Others/",
      gene,
      "_Plot.png"
    ),
    width = 7,
    height = 5
  )
}
########################################



# 1.3) Mitochondrial Ribosomal Proteins:

# 1.3.1) Small Sub unit (MRPS29 not present in counts data frame):
########################################
glist <- c(
  "MRPS2", "MRPS5", "MRPS6", "MRPS7", "MRPS9", "MRPS10", "MRPS11", "MRPS12", 
  "MRPS14", "MRPS15", "MRPS16", "MRPS17", "MRPS18A", "MRPS18B", "MRPS18C", 
  "MRPS22", "MRPS23", "MRPS24", "MRPS25", "MRPS26", "MRPS27", "MRPS28", "MRPS29", 
  "MRPS30", "MRPS31", "MRPS33", "MRPS34", "MRPS35"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}

glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)
table(res$Group)

for (gene in glist) {
  p_value <- Results[gene, ]$fdr
  
  if (is_scientific(p_value)) {
    annots <- signif(p_value, digits = 3)
  }
  else {
    annots <- round(p_value, 4)
  }
  
  ggplot(res, aes(x = AAA, y = res[, gene], fill = AAA)) +
    geom_boxplot() +
    ggtitle(gene) +
    ylab("Expression (Normalized TPM)") +
    labs(x = gene) +
    theme(
      plot.title = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.margin = margin(t = 0),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.title.x = element_text(
        size = 12,
        angle = 45,
        vjust = 0.5
      ),
      panel.background = element_blank()
    ) +
    scale_fill_manual(values = c("turquoise2", "salmon2"),
                      labels = c("Controls", "AAA")) +
    geom_signif(
      comparisons = list(c("Control", "Case")),
      map_signif_level = TRUE,
      textsize = 4,
      test = "t.test",
      annotations = annots
    )
  ggsave(
    paste0(
      "C://Users/Gerard/Desktop/Plots Hui/Results/Plots_Ribosomal_Small/",
      gene,
      "_Plot.png"
    ),
    width = 7,
    height = 5
  )
}
########################################


# 1.3.2) Large Sub unit:
########################################
glist <- c(
  "MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL9", "MRPL10", "MRPL11", "MRPL12", 
  "MRPL13", "MRPL14", "MRPL15", "MRPL16", "MRPL17", "MRPL18", "MRPL19", "MRPL20", 
  "MRPL21", "MRPL22", "MRPL23", "MRPL24", "MRPL27", "MRPL28", "MRPL30", "MRPL32", 
  "MRPL33", "MRPL34", "MRPL35", "MRPL36", "MRPL37", "MRPL38", "MRPL39", "MRPL40", 
  "MRPL41", "MRPL42", "MRPL43", "MRPL44", "MRPL45"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}

glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)
table(res$Group)

for (gene in glist) {
  p_value <- Results[gene, ]$fdr
  
  if (is_scientific(p_value)) {
    annots <- signif(p_value, digits = 3)
  }
  else {
    annots <- round(p_value, 4)
  }
  
  ggplot(res, aes(x = AAA, y = res[, gene], fill = AAA)) +
    geom_boxplot() +
    ggtitle(gene) +
    ylab("Expression (Normalized TPM)") +
    labs(x = gene) +
    theme(
      plot.title = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.margin = margin(t = 0),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.title.x = element_text(
        size = 12,
        angle = 45,
        vjust = 0.5
      ),
      panel.background = element_blank()
    ) +
    scale_fill_manual(values = c("turquoise2", "salmon2"),
                      labels = c("Controls", "AAA")) +
    geom_signif(
      comparisons = list(c("Control", "Case")),
      map_signif_level = TRUE,
      textsize = 4,
      test = "t.test",
      annotations = annots
    )
  ggsave(
    paste0(
      "C://Users/Gerard/Desktop/Plots Hui/Results/Plots_Ribosomal_Large/",
      gene,
      "_Plot.png"
    ),
    width = 7,
    height = 5
  )
}
########################################



# 1.4) Mitochondrial Transport Related Proteins:

# 1.4.1) Import (TOMM70A not present in counts dataframe) (TOMM40L is duplicated):
########################################
glist <- c(
  "TOMM20", "TOMM22", "TOMM40", "TOMM40L", "TOMM70A", "TOMM5", "TOMM6", "TOMM7", 
  "TIMM8A", "TIMM8B", "TIMM9", "TIMM10", "TIMM13", "TIMM17A", "TIMM17B", "TIMM22", 
  "TIMM23", "TIMM44", "SAMM50", "TOMM34", "DNAJC19", "PMPCA", "PMPCB", 
  "HSPA9", "DNAJA3", "HSPA8"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}

glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)
table(res$Group)

for (gene in glist) {
  p_value <- Results[gene, ]$fdr
  
  if (is_scientific(p_value)) {
    annots <- signif(p_value, digits = 3)
  }
  else {
    annots <- round(p_value, 4)
  }
  
  ggplot(res, aes(x = AAA, y = res[, gene], fill = AAA)) +
    geom_boxplot() +
    ggtitle(gene) +
    ylab("Expression (Normalized TPM)") +
    labs(x = gene) +
    theme(
      plot.title = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.margin = margin(t = 0),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.title.x = element_text(
        size = 12,
        angle = 45,
        vjust = 0.5
      ),
      panel.background = element_blank()
    ) +
    scale_fill_manual(values = c("turquoise2", "salmon2"),
                      labels = c("Controls", "AAA")) +
    geom_signif(
      comparisons = list(c("Control", "Case")),
      map_signif_level = TRUE,
      textsize = 4,
      test = "t.test",
      annotations = annots
    )
  ggsave(
    paste0(
      "C://Users/Gerard/Desktop/Plots Hui/Results/Plots_Transport_Import/",
      gene,
      "_Plot.png"
    ),
    width = 7,
    height = 5
  )
}
########################################


# 1.4.2) Export (SLC25A21, SLC25A31, SLC25A41 are removed for low expression):
########################################
glist <- c(
  "ABCB6", "ABCB7", "ABCB8", "ABCB10", "ABCC1", "ABCC9", "SLC25A1", "SLC25A3", 
  "SLC25A4", "SLC25A5", "SLC25A6", "SLC25A10", "SLC25A11", "SLC25A12", "SLC25A13", 
  "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17", "SLC25A19", "SLC25A20", "SLC25A21", 
  "SLC25A22", "SLC25A23", "SLC25A24", "SLC25A26", "SLC25A27", "SLC25A28", "SLC25A29", 
  "SLC25A30", "SLC25A31", "SLC25A32", "SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", 
  "SLC25A37", "SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42", "SLC25A43", 
  "SLC25A44", "SLC25A45", "SLC25A46", "TSPO", "VDAC1", "VDAC2", "VDAC3"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}

glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)
table(res$Group)

for (gene in glist) {
  p_value <- Results[gene, ]$fdr
  
  if (is_scientific(p_value)) {
    annots <- signif(p_value, digits = 3)
  }
  else {
    annots <- round(p_value, 4)
  }
  
  ggplot(res, aes(x = AAA, y = res[, gene], fill = AAA)) +
    geom_boxplot() +
    ggtitle(gene) +
    ylab("Expression (Normalized TPM)") +
    labs(x = gene) +
    theme(
      plot.title = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.margin = margin(t = 0),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.title.x = element_text(
        size = 12,
        angle = 45,
        vjust = 0.5
      ),
      panel.background = element_blank()
    ) +
    scale_fill_manual(values = c("turquoise2", "salmon2"),
                      labels = c("Controls", "AAA")) +
    geom_signif(
      comparisons = list(c("Control", "Case")),
      map_signif_level = TRUE,
      textsize = 4,
      test = "t.test",
      annotations = annots
    )
  ggsave(
    paste0(
      "C://Users/Gerard/Desktop/Plots Hui/Results/Plots_Transport_Export/",
      gene,
      "_Plot.png"
    ),
    width = 7,
    height = 5
  )
}
########################################





# 2) Mitophagy clustering heatmap: (Park2 changed to PARK2) (Prohibitn2 changed to PHB2)

# 2.1 Generalized mitophagy related genes :
########################################
glist <- c(
  "PARK2", "PINK1", "USP30", "BNIP3", "BNIP3L", "FUNDC1", "BECN1", "AMBRA1",
  "OPTN", "PHB2", "PPARGC1A", "ATG5", "ATG7", "MAP1LC3B"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}


glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)

gene_data <- res[, -ncol(res)] # All columns except the last one which is 'Group'

# Extract the group information
group_info <- as.factor(c(rep("AAA", 96), rep("Control", 44)))

# Create a color annotation for the groups
annotation <- data.frame(Group = factor(group_info))
rownames(annotation) <- rownames(res)

# Define custom colors for the heatmap:
ann_colors <- list(Group = c("AAA" = "salmon2", "Control" = "turquoise2"))

# Transpose the data so that genes are in rows and samples in columns
gene_data_t <- t(gene_data)

# Save the heatmap to a PNG file
png("C://Users/Gerard/Desktop/Plots Hui/Results/HeatMaps/Generalized_Mitophagy.png", width = 4000, height = 3000, res = 300)
pheatmap(
  gene_data_t,
  scale = "row",
  annotation_col = annotation,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 8
)
dev.off()
########################################



# 2.2 Specific mitophagy related genes :
########################################
glist <- c(
  "PARK2", "PINK1", "USP30", "BNIP3", "BNIP3L", "FUNDC1", "PHB2"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}


glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)

gene_data <- res[, -ncol(res)] # All columns except the last one which is 'Group'

# Extract the group information
group_info <- as.factor(c(rep("AAA", 96), rep("Control", 44)))

# Create a color annotation for the groups
annotation <- data.frame(Group = factor(group_info))
rownames(annotation) <- rownames(res)

# Define custom colors for the heatmap:
ann_colors <- list(Group = c("AAA" = "salmon2", "Control" = "turquoise2"))

# Transpose the data so that genes are in rows and samples in columns
gene_data_t <- t(gene_data)

# Save the heatmap to a PNG file
png("C://Users/Gerard/Desktop/Plots Hui/Results/HeatMaps/Specific_Mitophagy.png", width = 4000, height = 3000, res = 300)
pheatmap(
  gene_data_t,
  scale = "row",
  annotation_col = annotation,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 8
)
dev.off()
########################################



# 2.3 Specific mitophagy related genes, with potential upstream genes:
########################################
glist <- c(
  "PARK2", "PINK1", "USP30", "BNIP3", "BNIP3L", "FUNDC1", "PHB2", "BECN1", "AMBRA1", "OPTN"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}


glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)

gene_data <- res[, -ncol(res)] # All columns except the last one which is 'Group'

# Extract the group information
group_info <- as.factor(c(rep("AAA", 96), rep("Control", 44)))

# Create a color annotation for the groups
annotation <- data.frame(Group = factor(group_info))
rownames(annotation) <- rownames(res)

# Define custom colors for the heatmap:
ann_colors <- list(Group = c("AAA" = "salmon2", "Control" = "turquoise2"))

# Transpose the data so that genes are in rows and samples in columns
gene_data_t <- t(gene_data)

# Save the heatmap to a PNG file
png("C://Users/Gerard/Desktop/Plots Hui/Results/HeatMaps/Specific_Mitophagy_Upstream.png", width = 4000, height = 3000, res = 300)
pheatmap(
  gene_data_t,
  scale = "row",
  annotation_col = annotation,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 8
)
dev.off()
########################################



# 2.4 PINK1/PARK2 Pathway:
########################################
glist <- c(
  "PARK2", "PINK1", "USP30"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}


glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)

gene_data <- res[, -ncol(res)] # All columns except the last one which is 'Group'

# Extract the group information
group_info <- as.factor(c(rep("AAA", 96), rep("Control", 44)))

# Create a color annotation for the groups
annotation <- data.frame(Group = factor(group_info))
rownames(annotation) <- rownames(res)

# Define custom colors for the heatmap:
ann_colors <- list(Group = c("AAA" = "salmon2", "Control" = "turquoise2"))

# Transpose the data so that genes are in rows and samples in columns
gene_data_t <- t(gene_data)

# Save the heatmap to a PNG file
png("C://Users/Gerard/Desktop/Plots Hui/Results/HeatMaps/PINK1_PARK2_Pathway.png", width = 4000, height = 3000, res = 300)
pheatmap(
  gene_data_t,
  scale = "row",
  annotation_col = annotation,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 8
)
dev.off()
########################################



# 2.5 PINK1/PARK2 Pathway, with potential upstream genes:
########################################
glist <- c(
  "PARK2", "PINK1", "USP30", "BECN1", "AMBRA1", "OPTN"
)

# Check presence of gene:
for (gene in glist) {
  if (!(gene %in% rownames(counts))) {
    print(paste("Gene", gene, "is not present in the counts dataframe."))
    if (gene %in% rownames(counts.rm)) {
      print(paste("Gene", gene, "is removed for low expression."))
    }
  }
}


glist <- glist[glist %in% rownames(counts)]
glist


# Create empty data frame:
res <- data.frame(matrix(ncol = length(glist), nrow = 140))
colnames(res) <- glist
rownames(res) <- colnames(counts)
head(res)

res$Group <- AAA

# Fill the NAs of the res data frame:
for (gene in glist) {
  try(res_gene <- resid(lm(counts[gene,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)))
  res[[gene]] <- res_gene
}

head(res)

gene_data <- res[, -ncol(res)] # All columns except the last one which is 'Group'

# Extract the group information
group_info <- as.factor(c(rep("AAA", 96), rep("Control", 44)))

# Create a color annotation for the groups
annotation <- data.frame(Group = factor(group_info))
rownames(annotation) <- rownames(res)

# Define custom colors for the heatmap:
ann_colors <- list(Group = c("AAA" = "salmon2", "Control" = "turquoise2"))

# Transpose the data so that genes are in rows and samples in columns
gene_data_t <- t(gene_data)

# Save the heatmap to a PNG file
png("C://Users/Gerard/Desktop/Plots Hui/Results/HeatMaps/PINK1_PARK2_Pathway_Upstream.png", width = 4000, height = 3000, res = 300)
pheatmap(
  gene_data_t,
  scale = "row",
  annotation_col = annotation,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 8
)
dev.off()
########################################
