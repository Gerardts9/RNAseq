library(openxlsx);library(dplyr)

Results_CC <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/AllRes/Cases_Controls.xlsx")
Results_DIA <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/AllRes/Diameter.xlsx")
Results_SYM <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/AllRes/Symptoms.xlsx")

sig.genes <- rownames(Results_CC[Results_CC$fdr < 0.0000000005,])

common <- intersect(intersect(rownames(Results_CC), rownames(Results_DIA)), rownames(Results_SYM))

Results_CC_2 <- Results_CC[rownames(Results_CC) %in% common, ]
Results_DIA_2 <- Results_DIA[rownames(Results_DIA) %in% common, ]
Results_SYM_2 <- Results_SYM[rownames(Results_SYM) %in% common, ]

hist(Results_CC_2$Beta[Results_CC_2$Beta > -20 & Results_CC_2$Beta < 20])
hist(Results_DIA_2$Beta[Results_DIA_2$Beta > -0.7 & Results_DIA_2$Beta < 0.7])
hist(Results_SYM_2$Beta[Results_SYM_2$Beta > -10 & Results_SYM_2$Beta < 10])



# P-values*Beta Plots:
############################
Results_CC_2 <- Results_CC_2 %>%
  mutate(plot = ifelse(Beta > 0, -log10(pvalue), log10(pvalue)))

Results_DIA_2 <- Results_DIA_2 %>%
  mutate(plot = ifelse(Beta > 0, -log10(pvalue), log10(pvalue)))

Results_SYM_2 <- Results_SYM_2 %>%
  mutate(plot = ifelse(Beta > 0, -log10(pvalue), log10(pvalue)))

res <- data.frame(
  CC = Results_CC_2$pvalue,
  DIA = Results_DIA_2$pvalue,
  SYM = Results_SYM_2$pvalue,
  row.names = common)

res.beta <- data.frame(
  CC = Results_CC_2$Beta,
  DIA = Results_DIA_2$Beta,
  SYM = Results_SYM_2$Beta,
  row.names = common)

mean(res$CC) # -1.82
mean(res$DIA) # 0.26
mean(res$SYM) # -0.37

range(res$CC) # -46.38  34.62
range(res$DIA) # -13.31  12.76
range(res$SYM) # -10.66  17.28

plot(-log10(res$CC), -log10(res$DIA), xlim = c(0,20), ylim = c(0,20))
plot(res.beta$CC, res.beta$DIA, xlim = c(-20, 20), ylim = c(-20, 20))


plot(-log10(res$CC), -log10(res$SYM), xlim = c(0,20), ylim = c(0,20))
plot(res.beta$CC, res.beta$SYM, xlim = c(-100, 100), ylim = c(-100, 100))
############################



# Top Case Control:
############################
cc.sig <- rownames(Results_CC[Results_CC$fdr < 0.0000000000005, ])

head(Results_CC)
head(Results_DIA)
head(Results_SYM)

Results_CC_2 <- Results_CC[rownames(Results_CC) %in% cc.sig, ]
Results_DIA_2 <- Results_DIA[rownames(Results_DIA) %in% cc.sig, ]
Results_SYM_2 <- Results_SYM[rownames(Results_SYM) %in% cc.sig, ]

combined_df <- merge(Results_CC_2, Results_DIA_2, by = 0 )
combined_df <- merge(combined_df, Results_SYM, by.x = "Row.names", by.y = 0)
combined_df <- combined_df %>% remove_rownames %>% column_to_rownames(var="Row.names")

head(combined_df)
names(combined_df)[names(combined_df) == c("Beta.x")] <- "Case"
names(combined_df)[names(combined_df) == c("Beta.y")] <- "Diameter"
names(combined_df)[names(combined_df) == c("Beta")] <- "Symptomatology"

png("C://Users/Gerard/Desktop/AAA/RNAseq/Plots/Heatmap/cc.sig.png", width = 750, height = 750)
pheatmap(combined_df[,c("Case","Diameter","Symptomatology")], main = "Heatmap of Differential Expression Results")
dev.off()
############################



# Diameter Significant:
############################
dia.sig <- rownames(Results_DIA[Results_DIA$fdr < 0.05, ])

head(Results_CC)
head(Results_DIA)
head(Results_SYM)

Results_CC_2 <- Results_CC[rownames(Results_CC) %in% dia.sig, ]
Results_DIA_2 <- Results_DIA[rownames(Results_DIA) %in% dia.sig, ]
Results_SYM_2 <- Results_SYM[rownames(Results_SYM) %in% dia.sig, ]

combined_df <- merge(Results_CC_2, Results_DIA_2, by = 0 )
combined_df <- merge(combined_df, Results_SYM, by.x = "Row.names", by.y = 0)
combined_df <- combined_df %>% remove_rownames %>% column_to_rownames(var="Row.names")

head(combined_df)
names(combined_df)[names(combined_df) == c("Beta.x")] <- "Case"
names(combined_df)[names(combined_df) == c("Beta.y")] <- "Diameter"
names(combined_df)[names(combined_df) == c("Beta")] <- "Symptomatology"

png("C://Users/Gerard/Desktop/AAA/RNAseq/Plots/Heatmap/Dia.sig.png", width = 750, height = 750)
pheatmap(combined_df[,c("Case","Diameter","Symptomatology")], main = "Heatmap of Differential Expression Results")
dev.off()
############################




# Symptom Significant:
############################
sym.sig <- rownames(Results_SYM[Results_SYM$fdr < 0.05, ])

head(Results_CC)
head(Results_DIA)
head(Results_SYM)

Results_CC_2 <- Results_CC[rownames(Results_CC) %in% sym.sig, ]
Results_DIA_2 <- Results_DIA[rownames(Results_DIA) %in% sym.sig, ]
Results_SYM_2 <- Results_SYM[rownames(Results_SYM) %in% sym.sig, ]

combined_df <- merge(Results_CC_2, Results_DIA_2, by = 0 )
combined_df <- merge(combined_df, Results_SYM, by.x = "Row.names", by.y = 0)
combined_df <- combined_df %>% remove_rownames %>% column_to_rownames(var="Row.names")

head(combined_df)
names(combined_df)[names(combined_df) == c("Beta.x")] <- "Case"
names(combined_df)[names(combined_df) == c("Beta.y")] <- "Diameter"
names(combined_df)[names(combined_df) == c("Beta")] <- "Symptomatology"

png("C://Users/Gerard/Desktop/AAA/RNAseq/Plots/Heatmap/Sym.sig.png", width = 750, height = 750)
pheatmap(combined_df[,c("Case","Diameter","Symptomatology")], main = "Heatmap of Differential Expression Results")
dev.off()
############################










