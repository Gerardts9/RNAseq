library(data.table);library(edgeR);library(tidyverse);library(openxlsx);library(preprocessCore);library("org.Hs.eg.db");library(stats);library(qvalue);library(ggpubr);library(VennDiagram)
library(readxl);library(tidyr)

lm22 <- data.frame(read_xls("C://Users/Gerard/Downloads/41592_2015_BFnmeth3337_MOESM207_ESM.xls", sheet = 1, skip = 13))

head(lm22)

#fwrite(lm22[,3:ncol(lm22)], "C://Users/Gerard/Desktop/AAA/RNAseq/Cibersort/lm22.txt", sep = "\t")


ann <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/Annotation/gencode.v26.annotation.fixed.gtf.gz")
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

# Example for 1 row:
res <- resid(lm(tpm_matrix[1,] ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO))


# Define a function to calculate residuals for a row:
calculate_residuals <- function(row) {
  model <- lm(row ~ tec$FC + tec$Lane + tec$GC_Mean + tec$RIN + tec$DV200 + tec$Qubit + bio$age + bio$SEXO)
  residuals <- resid(model)
  return(residuals)
}

# Apply the function to each row of the tpm_matrix using apply():
residuals <- data.frame(t(apply(counts.tpm, 1, calculate_residuals)))

residuals[1:5,1:5]
res[1:5]

residuals2 <- merge(residuals, ann[,c("gene_id","gene_name")], by.x = 0, by.y = "gene_id")
names(residuals2)[names(residuals2) == "Row.names"] <- "gene_id"
residuals2 <- residuals2[!duplicated(residuals2$gene_name),]
residuals2 <- residuals2 %>% remove_rownames %>% column_to_rownames(var="gene_name")
residuals2 <- residuals2[,-1]

head(residuals2)

residuals2[1:5,1:5]
res[1:5]

#fwrite(residuals2, "C://Users/Gerard/Desktop/AAA/RNAseq/Cibersort/resid_counts.txt", sep = "\t", row.names = T)



# Analyse results:
res <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/Cibersort/CIBERSORTx_Job5_Adjusted.txt")
res <- res %>% remove_rownames %>% column_to_rownames(var="Mixture")
res <- res[,1:22]
head(res)

ncol(res)

rowSums(res)

df_cases <- res[1:96,]
df_controls <- res[97:140,]

head(df_cases)
head(df_controls)

rowSums(df_cases)
rowSums(df_controls)

column_mean_cases <- 100*apply(df_cases, 2, mean)
sum(column_mean_cases)
sorted_column_mean_cases <- sort(column_mean_cases, decreasing = TRUE)
top_three_cell_types_cases <- head(sorted_column_mean_cases, n = 3)
less_three_cell_types_cases <- tail(sorted_column_mean_cases, n = 3)

sum(sorted_column_mean_cases)

column_mean_controls <- 100*apply(df_controls, 2, mean)
sorted_column_mean_controls <- sort(column_mean_controls, decreasing = TRUE)
top_three_cell_types_controls <- head(sorted_column_mean_controls, n = 3)
less_three_cell_types_controls <- tail(sorted_column_mean_controls, n = 3)

sum(sorted_column_mean_controls)

res <- res %>%
  mutate(group = ifelse(row_number() <= 96, "Cases", "Controls"))


library(reshape2)
df.m <- melt(res, id.vars = "group")
head(df.m)

class(df.m$value)

ggplot(df.m, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free") +
  labs(x = "group", y = "Expression Level") +
  theme_bw() +
  theme(legend.position = "none")

t.test_results <- df.m %>%
  group_by(variable) %>%
  summarize(p_value = t.test(value ~ group)$p.value)

t.test_results[t.test_results$p_value < 0.05/22,]

t.test_results$p_value <- as.character(round(t.test_results$p_value, 4))

t.test_results[t.test_results$variable %in% "Dendritic cells activated",]$p_value <- "3.17e-5"
t.test_results[t.test_results$variable %in% "T cells CD8",]$p_value <- "3.83e-4"

t.test_results$p_value

# CD8 T-cells. HIGHER IN CASES.
# NK resting cells. HIGHER IN CONTROLS.
# Dendritic activated cells. HIGHER IN CONTROLS.

top_three_cell_types_cases
less_three_cell_types_cases

top_three_cell_types_controls
less_three_cell_types_controls

top_three_cell_types

ggplot(df.m, aes(x = group, y = value)) +
  geom_boxplot() +
  labs(x = "group", y = "Expression Level") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~ variable, scales = "free") +
  geom_text(
    data = t.test_results,
    aes(x = 1, y = max(df.m$value), 
        label = paste("p-value =", p_value)),
    hjust = 0, vjust = 1
  )


df.m <- left_join(df.m, t.test_results, by = "variable")


ggplot(df.m, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  labs(x = "group", y = "Expression Level") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~ variable, scales = "free") +
  geom_text(
    data = df.m,
    aes(x = 1, y = max(value), 
        label = paste("p-value =", p_value)),
    hjust = 0, vjust = 1
  )


ggsave("C://Users/Gerard/Desktop/AAA/RNAseq/Cibersort/Plot.png", width = 16, height = 14)













mpg_plot <- ggplot(mpg, aes(x = displ, y = hwy, fill = manufacturer )) +
  geom_point() +
  facet_grid(. ~ drv)

# A data frame with labels for each facet
f_labels <- data.frame(drv = c("4", "f", "r"), label = c("4wd", "Front", "Rear"))

mpg_plot +
  geom_text(x = 6, y = 40, aes(label = label), data = f_labels)                                 
                                    
                                    
                                    
                  