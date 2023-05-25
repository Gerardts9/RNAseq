library(data.table);library(dplyr);library(tidyverse);library(openxlsx);library(ggpubr);library(preprocessCore);library(readxl);library(qvalue);library(dplyr)


ann <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/Annotation/gencode.v26.annotation.fixed.gtf.gz")
ann$gene_id <- substr(ann$gene_id, 1, 15)

# Load Counts:
########################################
counts.all <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/Gene_counts/Counts.RSEM.All.txt")
counts.all <- merge(counts.all, ann[,c("gene_id","gene_name")], by.x = 0, by.y = "gene_id")
counts.all <- counts.all[,-1]
counts.all <- counts.all[!duplicated(counts.all$gene_name),]
counts.all <- counts.all %>% remove_rownames %>% column_to_rownames(var="gene_name")
counts.all <- as.matrix(normalize.quantiles(as.matrix(counts.all), keep.names = T))

counts.cases <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/Gene_counts/Counts.RSEM.Cases.txt")
counts.cases <- merge(counts.cases, ann[,c("gene_id","gene_name")], by.x = 0, by.y = "gene_id")
counts.cases <- counts.cases[,-1]
counts.cases <- counts.cases[!duplicated(counts.cases$gene_name),]
counts.cases <- counts.cases %>% remove_rownames %>% column_to_rownames(var="gene_name")
counts.cases <- as.matrix(normalize.quantiles(as.matrix(counts.cases), keep.names = T))
########################################


# Load biological/technical data:
########################################
tec2 <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/Tables/Technical_Complete.xlsx")
tec2 <- tec2[tec2$Seq_ID %in% colnames(counts.all),]
bio2 <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/20221115_Datos Pacientes.xlsx")
bio2 <- bio2[bio2$Muestra %in% colnames(counts.all),]

tec1 <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/Tables/First Analysis/Technical_CompleteV2.xlsx")
tec1 <- tec1[tec1$SampleID %in% colnames(counts.cases),]
bio1 <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/20221115_Datos Pacientes.xlsx")
bio1 <- bio1[bio1$Muestra %in% colnames(counts.cases),]
########################################


# Load results:
########################################
cc <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Cases_Controls.txt")
dia <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Diameter.txt")
sym <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Symptoms.txt")

alt <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/SUPPA/AAA_Sig_Results.txt")
alt$event <- rownames(alt)
alt$gene_id2 <- substr(rownames(alt), 1, 15)
alt <- merge(alt, ann, by = "gene_id2")

med <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 1)
med$Gene.Symbol <- gsub("-", "", med$Gene.Symbol)

adv <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 2)
adv$Gene.Symbol <- gsub("-", "", adv$Gene.Symbol)
########################################


intersect(cc$V1, intersect(dia$V1, sym$V1)) # LAMA2

intersect(cc$V1, intersect(dia$V1, intersect(sym$V1, alt$gene_name)))

intersect(cc$V1, alt$gene_name)


# LAMA2:
which(cc$V1 == "LAMA2") # 3507
which(dia$V1 == "LAMA2") # 29
which(sym$V1 == "LAMA2") # 1

lama <- ann[ann$gene_name == "LAMA2",]$gene_id


# By Case-Control:
########################################
cc[cc$V1 == "LAMA2",]

mean(counts.all["LAMA2", 1:96])
mean(counts.all["LAMA2", 97:140])

res1 <- resid(lm(counts.all["LAMA2",] ~ tec2$FC + tec2$Lane + tec2$GC_Mean + tec2$RIN + tec2$DV200 + tec2$Qubit + bio2$age + bio2$SEXO))

res2 <- lm(counts.all["LAMA2",] ~ tec2$FC + tec2$Lane + tec2$GC_Mean + tec2$RIN + tec2$DV200 + tec2$Qubit + bio2$age + bio2$SEXO)$residuals

identical(res1, res2)

res <-  data.frame(Expression = res1)

ggplot(res, aes(x = bio2$Status, y = Expression, fill = bio2$Status)) + geom_boxplot(outlier.shape = NA) + ggtitle("LAMA2") +
  ylab("Expression (Normalized TPM)") + xlab("Status") +
  geom_point(position = position_jitter(width = 0.2), color = "black", size = 2, alpha = 0.5) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  ) + scale_fill_discrete(name = "Status")

#ggsave("C://Users/Gerard/Desktop/AAA/RNAseq/Plots/Results/LAMA2_CaseControl.png", width = 7.5, height = 7.5)
########################################


# By diameter:
########################################
dia <- bio1[!is.na(bio1$aortc_diameter_mm),][,c("Muestra","aortc_diameter_mm")]

counts.cases2 <- counts.cases[,colnames(counts.cases) %in% dia$Muestra]
tec.dia <- tec1[tec1$SampleID %in% dia$Muestra,]
bio.dia <- bio1[bio1$Muestra %in% dia$Muestra,]

dim(counts.cases2)
dim(tec.dia)
dim(bio.dia)

dia <- dia %>% remove_rownames %>% column_to_rownames(var="Muestra")

res <- resid(lm(counts.cases2["LAMA2",] ~ tec.dia$Date + tec.dia$Batch + tec.dia$GC_Mean + tec.dia$RIN + 
                  tec.dia$DV200 + tec.dia$Qubit + bio.dia$age))

res <- data.frame(Expression = res,
                  Diameter = dia)

ggplot(res, aes(x = aortc_diameter_mm, y = Expression)) + geom_point(color = "orangered2") + ggtitle("LAMA2") + ylab ("Expression (Normalized TPM)") + xlab("Aortic Diameter (mm)") +
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
    axis.line = element_line(colour = "black")
  )

#ggsave(filename = "C://Users/Gerard/Desktop/AAA/RNAseq/Plots/Results/LAMA2_Diameter.png", width = 7.5, height = 7.5)
########################################


# By Symptoms:
########################################
no <- bio1[bio1$symptomatic == 0,]$Muestra
yes <- bio1[bio1$symptomatic == 1,]$Muestra

mean(counts.all["LAMA2", colnames(counts.all) %in% no])
mean(counts.all["LAMA2", colnames(counts.all) %in% yes])

sym <- bio1[!is.na(bio1$symptomatic),][,c("Muestra","symptomatic")]

counts.cases2 <- counts.cases[,colnames(counts.cases) %in% sym$Muestra]
tec.sym <- tec1[tec1$SampleID %in% sym$Muestra,]
bio.sym <- bio1[bio1$Muestra %in% sym$Muestra,]

dim(counts.cases2)
dim(tec.sym)
dim(bio.sym)

sym <- sym %>% remove_rownames %>% column_to_rownames(var="Muestra")

res <- resid(lm(counts.cases2["LAMA2",] ~ tec.sym$Date + tec.sym$Batch + tec.sym$GC_Mean + tec.sym$RIN + 
                  tec.sym$DV200 + tec.sym$Qubit + bio.sym$age))

res <- data.frame(Expression = res,
                  Type = sym)

res$symptomatic <- ifelse(res$symptomatic == 0, "Asymptomatic", "Symptomatic")
table(res$symptomatic)

ggplot(res, aes(x = res$symptomatic, y = Expression, fill = res$symptomatic)) + geom_boxplot(outlier.shape = NA) + ggtitle("LAMA2") +
  ylab("Expression (Normalized TPM)") + xlab("Symptoms") +
  geom_point(position = position_jitter(width = 0.2), color = "black", size = 2, alpha = 0.5) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  ) + scale_fill_discrete(name = "Type")

#ggsave("C://Users/Gerard/Desktop/AAA/RNAseq/Plots/Results/LAMA2_Symptoms.png", width = 7.5, height = 7.5)
########################################


# Run 1.3.Manual_Analysis_RSEM_Cases.Rmd diameter results:

hist(Results[rownames(Results) %in% med$Gene.Symbol,]$pvalue)
hist(Results[rownames(Results) %in% adv$Gene.Symbol,]$pvalue)


# pi1 Media in Diameter results:
1-qvalue(p = Results[rownames(Results) %in% med$Gene.Symbol,]$pvalue)$pi0 # 0.37

# pi1 Adventitia in Diameter results:
1-qvalue(p = Results[rownames(Results) %in% adv$Gene.Symbol,]$pvalue)$pi0 # 0.26

# pi CasesControls in Diameter results:
1-qvalue(p = Results[rownames(Results) %in% cc$V1,]$pvalue)$pi0 # 0.39

# 37% of the DEG between Media cases and controls show some sort of significance in DEG by diameter.

# 26% of the DEG between Adventitia cases and controls show some sort of significance in DEG by diameter.

# 40% of the DEG between cases and controls show some sort of significance in DEG by diameter.



# Contingency table CaseControl vs Splicing:
res <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/SUPPA/AAA_Results.dpsi")
res2 <- na.omit(res)
res2$gene_id <- substr(rownames(res2), 1, 15)
res2 <- merge(res2, ann, by = "gene_id")

sig.alt <- unique(res2[res2$AAA_Cases.AAA_Controls_dPSI > 0.1 & res2$AAA_Cases.AAA_Controls_p.val < 0.05, ]$gene_id)
no.sig.alt <- unique(res2[!res2$gene_id %in% sig.alt,]$gene_id)


Results$gene_id <- substr(Results$gene_id, 1, 15)

sig.cc <- unique(Results[Results$fdr < 0.05,]$gene_id)
no.sig.cc <- unique(Results[!Results$gene_id %in% sig.cc,]$gene_id)

length(unique(res2$gene_id)) # 9278
length(sig.alt) + length(no.sig.alt)

length(unique(Results$gene_id)) # 10573
length(sig.cc) + length(no.sig.cc)


length(intersect(sig.cc, sig.alt))
length(intersect(sig.cc, no.sig.alt))
length(intersect(sig.alt, no.sig.cc))
length(intersect(no.sig.alt, no.sig.cc))

