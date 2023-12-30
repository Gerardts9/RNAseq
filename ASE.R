# Create plot to represent the imbalance:
######################
ff <-  list.files("/home/gerard/AAA/ASE/Results/", pattern = "*_case_gene_ae.txt", recursive = T, full.names = T)

for (i in 1:length(ff)) {
  folder <- gsub("_case_gene_ae.txt", "", basename(ff[i]))
  name <- basename(ff[i])
  
  # load haplotypic counts
  phaser_test_case_gene_ae = read.delim(paste0("/home/gerard/AAA/ASE/Results/", folder, "/", name))
  # select only genes with sufficient coverage
  cov8 = subset(phaser_test_case_gene_ae,
                phaser_test_case_gene_ae$totalCount >= 8)
  # perform a binomial test for deviation from 0.5
  cov8$binom_p = apply(cov8[, c("aCount", "bCount")], 1, function(x)
    binom.test(x[1], x[1] + x[2], p = 0.5)$p.value)
  # perform multiple testing correction with FDR
  cov8$binom_q = p.adjust(cov8$binom_p, method = "fdr")
  # plot haplotype A versus haplotype, highlight sites with significant imbalance (FDR < 5%)
  
  png(paste0(
    "/home/gerard/AAA/ASE/Results/",
    folder,
    "/",
    folder,
    "_Results.png"
  ))
  plot(
    cov8$bCount,
    cov8$aCount,
    log = "xy",
    col = (cov8$binom_q < 0.05) + 1,
    ylab = "Haplotype B Count",
    xlab = "Haplotype A Count"
  )
  abline(0, 1, col = "grey")
  legend(
    "topleft",
    c("No significant imbalance", "Significant imbalance"),
    pch = c(15, 15),
    col = c(1, 2)
  )
  dev.off()
}
######################


# QQ plot of p-values:
######################
ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}

for (i in 1:length(ff)) {
  folder <- gsub("_case_gene_ae.txt", "", basename(ff[i]))
  name <- basename(ff[i])

  # load haplotypic counts
  phaser_test_case_gene_ae = read.delim(paste0("/home/gerard/AAA/ASE/Results/", folder, "/", name))
  # select only genes with sufficient coverage
  cov8 = subset(phaser_test_case_gene_ae,
                phaser_test_case_gene_ae$totalCount >= 8)
  # perform a binomial test for deviation from 0.5
  cov8$binom_p = apply(cov8[, c("aCount", "bCount")], 1, function(x)
    binom.test(x[1], x[1] + x[2], p = 0.5)$p.value)
  # perform multiple testing correction with FDR
  cov8$binom_q = p.adjust(cov8$binom_p, method = "fdr")
  # plot haplotype A versus haplotype, highlight sites with significant imbalance (FDR < 5%)
  
png(paste0("/home/gerard/AAA/ASE/Results/", folder,"/", folder, "_QQPlot.png"))
ggd.qqplot(cov8[cov8$binom_p>0,]$binom_p)
dev.off()
}
######################


library(data.table);library(tidyr);library(tidyverse);library(ggplot2);library(ggbiplot);library(openxlsx);
library(dplyr);library(biomaRt);library(clusterProfiler);library(aPEAR);library(org.Hs.eg.db)


# Function to split values and create new columns for each gene
split_alleles <- function(df) {
  genes <- names(df)
  for (gene in genes) {
    alleles <- strsplit(as.character(df[[gene]]), "|", fixed = TRUE)
    a_counts <- sapply(alleles, function(x) as.numeric(strsplit(x, "|", fixed = TRUE)[[1]]))
    b_counts <- sapply(alleles, function(x) as.numeric(strsplit(x, "|", fixed = TRUE)[[2]]))
    
    colnames(df)[colnames(df) == gene] <- paste0(gene, "_aCount")
    df[[paste0(gene, "_aCount")]] <- a_counts
    
    df[[paste0(gene, "_bCount")]] <- b_counts
  }
  return(df)
}


# Allele 'a' and Allele 'b' counts for each control sample and gene:
######################
controls <- fread("/home/gerard/AAA/ASE/phASER_WASP_GTEx_v8_matrix.txt.gz")

# Select Artery Aorta samples: 
ss <- fread("/home/gerard/AAA/ASE/annotations_v8_metadata-files_GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
art.aor.ids <- ss[ss$SMTSD %in% "Artery - Aorta",]$SAMPID


selected_columns <- names(controls)[names(controls) %in% art.aor.ids]
selected_columns <- c("name", selected_columns)
controls <- controls[, selected_columns, with = FALSE]
controls <- controls %>% remove_rownames %>% column_to_rownames(var="name")
controls <- as.data.frame(t(controls))

controls_split <- split_alleles(controls)

controls_split2 <- controls_split[order(names(controls_split))]

# Check:
head(controls_split2[grep("ENSG00000227232.5", names(gtex2_split2))])

fwrite(controls_split2, "/home/gerard/AAA/ASE/Controls_Counts.txt", sep = "\t", row.names = TRUE)
######################



# Allele 'a' and Allele 'b' counts for each case sample and gene:
######################
ff <-  list.files("/home/gerard/AAA/ASE/Results/", pattern = "*_case_gene_ae.txt", recursive = T, full.names = T)

data_list <- lapply(ff, fread)
cases_combined <- bind_rows(data_list)


cases_combined <- cases_combined[,c("name", "bam", "aCount", "bCount")]
names(cases_combined) <- c("gene","sample","aCount","bCount")



cases_combined$combined <- paste(cases_combined$aCount, cases_combined$bCount, sep = "|")

cases_wide <- dcast(cases_combined, sample ~ gene, value.var = "combined")
cases_wide <- cases_wide %>% remove_rownames %>% column_to_rownames(var="sample")

cases_split <- split_alleles(cases_wide)

cases_split2 <- cases_split[order(names(cases_split))]


fwrite(cases_split2, "/home/gerard/AAA/ASE/Cases_Counts.txt", sep = "\t", row.names = TRUE)
######################



# Clustering with GTEx controls data:
######################
cases <- fread("/home/gerard/AAA/ASE/Cases_Counts.txt")
cases <- cases %>% remove_rownames %>% column_to_rownames(var="V1")

controls <- fread("/home/gerard/AAA/ASE/Controls_Counts.txt")
controls <- controls %>% remove_rownames %>% column_to_rownames(var="V1")

common <- intersect(names(cases),names(controls))

cases <- cases[common]
controls <- controls[common]

if (identical(names(cases),names(controls)) != TRUE){warning("Genes are not the same")}

all <- rbind(cases,controls)
all[1:5,1:5]

aaa <- as.factor(c(rep(1,nrow(cases)),rep(0, nrow(controls))))

my_pca <- prcomp(all, center = T)
summaryPca <- summary(my_pca)


png("/home/gerard/AAA/ASE/Plots/cumsum.png")
plot(cumsum(my_pca$sdev^2 / sum(my_pca$sdev^2))[1:10], type="b", ylab = "Cummulative proportion of Variance Explained")
dev.off()

png(file = "/home/gerard/AAA/ASE/Plots/Scree_plot.png")
fviz_eig(my_pca)
dev.off()


png(file = "/home/gerard/AAA/ASE/Plots/Individuals_plot.png", width = 1500, height = 1500, res = 200)
fviz_pca_ind(
  my_pca,
  habillage = aaa,
  palette = c("#00AFBB", "#FC4E07"),
  repel = TRUE,
  addEllipses = TRUE,
  # Add ellipses
  ellipse.level = 0.95,
  geom="point"
)
dev.off()



summed_data <- combined_data %>%
  group_by(name) %>%
  summarise(across(c(aCount, bCount), sum, na.rm = TRUE))


dim(summed_data) # 58219    3

summed_data <- summed_data %>% remove_rownames %>% column_to_rownames(var="name")

fwrite(summed_data, "/home/gerard/AAA/ASE/Cases_Counts.txt", sep = "\t", row.names = TRUE)

summed_data <- fread("/home/gerard/AAA/ASE/Cases_Counts.txt")

######################



# Compare ratios between cases and controls:
######################
cases <- fread("/home/gerard/AAA/ASE/Cases_Counts.txt")
cases <- cases %>% remove_rownames %>% column_to_rownames(var="V1")
cases <- cases[,1:6]

controls <- fread("/home/gerard/AAA/ASE/Controls_Counts.txt")
controls <- controls %>% remove_rownames %>% column_to_rownames(var="V1")
controls <- controls[,1:6]

common <- intersect(names(cases),names(controls))

cases <- cases[common]
controls <- controls[common]

if (identical(names(cases),names(controls)) != TRUE){warning("Genes are not the same")}


# Find columns containing "_aCount":
aCount_cols <- grep("_aCount", colnames(cases))

# Extract gene names from column names:
gene_names <- sub("_aCount", "", colnames(cases)[aCount_cols])

# Loop through each gene and calculate the sum of counts for cases:
for (gene in gene_names) {
  # Select columns for aCount and bCount for the current gene:
  aCount_col <- paste(gene, "_aCount", sep = "")
  bCount_col <- paste(gene, "_bCount", sep = "")
  
  # Calculate the sum for each sample:
  total_counts <- cases[[aCount_col]] + cases[[bCount_col]]
  
  # Create a new column with the total counts for the current gene:
  cases[[paste0(gene, "_totalCount")]] <- total_counts
}

cases <- cases[,order(names(cases))]
cases <- cases[, -grep("bCount", colnames(cases))]


for (gene in gene_names) {
  # Select columns for aCount and totalCount for the current gene
  aCount_col <- paste(gene, "_aCount", sep = "")
  totalCount_col <- paste(gene, "_totalCount", sep = "")
  
  # Calculate the ratio for each sample
  cases[[paste0(gene, "_ratio")]] <-
    cases[[aCount_col]] / cases[[totalCount_col]]
}

cases[is.na(cases)] <- 0
cases <- data.frame(round(cases, 4))
cases <- cases[, order(names(cases))]
cases <- cases[, -grep("aCount", colnames(cases))]
cases <- cases[, -grep("totalCount", colnames(cases))]
colnames(cases) <- gsub("_ratio","",colnames(cases))

fwrite(cases, "/home/gerard/AAA/ASE/Cases_ratio.txt", sep = "\t")




# Loop through each gene and calculate the sum of counts for controls:
for (gene in gene_names) {
  # Select columns for aCount and bCount for the current gene
  aCount_col <- paste(gene, "_aCount", sep = "")
  bCount_col <- paste(gene, "_bCount", sep = "")
  
  # Calculate the sum for each sample
  total_counts <- controls[[aCount_col]] + controls[[bCount_col]]
  
  # Create a new column with the total counts for the current gene
  controls[[paste0(gene, "_totalCount")]] <- total_counts
}

controls <- controls[,order(names(controls))]
controls <- controls[, -grep("bCount", colnames(controls))]

for (gene in gene_names) {
  # Select columns for aCount and totalCount for the current gene
  aCount_col <- paste(gene, "_aCount", sep = "")
  totalCount_col <- paste(gene, "_totalCount", sep = "")
  
  # Calculate the ratio for each sample
  controls[[paste0(gene, "_ratio")]] <-
    controls[[aCount_col]] / controls[[totalCount_col]]
}

controls[is.na(controls)] <- 0
controls <- data.frame(round(controls, 4))
controls <- controls[, order(names(controls))]
controls <- controls[, -grep("aCount", colnames(controls))]
controls <- controls[, -grep("totalCount", colnames(controls))]
colnames(controls) <- gsub("_ratio","",colnames(controls))

#fwrite(controls, "/home/gerard/AAA/ASE/Controls_ratio.txt", sep = "\t")


cases <- fread("/home/gerard/AAA/ASE/Cases_ratio.txt")
controls <- fread("/home/gerard/AAA/ASE/Controls_ratio.txt")

# Get the list of column names (genes)
gene_names <- colnames(cases)

# Initialize an empty list to store test results
wilcoxon_results <- list()

# Loop through each gene and perform Wilcoxon test
for (gene in gene_names) {
  # Extract data for the current gene for cases and controls
  data_cases <- cases[[gene]]
  data_controls <- controls[[gene]]
  
  # Perform Wilcoxon test
  wilcox_result <- wilcox.test(data_cases, data_controls)$p.value
  
  # Store the test result in the list
  wilcoxon_results[[gene]] <- wilcox_result
}


p_values <- unlist(wilcoxon_results)
p_values[is.na(p_values)] <- 1

adjusted_p_values <- p.adjust(p_values, method = "fdr")

significant_genes <- gene_names[adjusted_p_values < 0.05]

length(significant_genes) # 1815


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_ids <- substr(significant_genes, 1, 15)

gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                   filters = "ensembl_gene_id", 
                   values = ensembl_ids, 
                   mart = ensembl)

gene_info[gene_info$external_gene_name == "",]
gene_info <- gene_info[gene_info$external_gene_name != "",]


# Enrichment analysis with significant genes:
gene_list <- gene_info$external_gene_name

entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO enrichment analysis:
go_enrich <- enrichGO(gene          = entrez_ids$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL", # Choose from "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable = T)

go_enrich

clusters <- findPathClusters(go_enrich@result, cluster = 'hier', minClusterSize = 15)
clusters$clusters
ase.clusters <- unique(clusters$clusters$Cluster)
######################










# Concatenate all files:
aes <- list.files("C://Users/Gerard/Desktop/AAA/RNAseq/ASE/", pattern = "case_gene_ae.txt", full.names = T)
over.8 <- c()

for (i in 1:length(aes)) {
  ff <- fread(aes[i])
  ff$Table <- i
  ff2 <- subset(ff, ff$totalCount >= 8)
  over.8 <- rbind(over.8, ff2)
}

# Perform a binomial test for deviation from 0.5:
over.8$binom_p <- apply(over.8[,c("aCount","bCount")], 1, function(x) binom.test(x[1],x[1]+x[2],p=0.5)$p.value)

# Perform multiple testing correction with FDR for the 12 individuals:
over.8$binom_q <- p.adjust(over.8$binom_p, method = "fdr")

groups <- split(over.8, over.8$Table)

sig.rows <- function(df, columna) {
  sum(df[[columna]] < 0.05, na.rm = TRUE)
}

sig.rows.table <- sapply(groups, function(df) {
  sig.rows(df, "binom_q")
})

mean(sig.rows.table)



# In average, there are 530 genes with allele specific expression in the 12 individuals with genotype information.



# Study overlap among significant genes between individuals:

sig.list <- lapply(groups, subset, binom_q < 0.05)
sig.list[[12]]


column_name <- "name"

percentage_overlaps <- c()
overlaps <- c()

for (i in 1:(length(sig.list)-1)) {
  for (j in (i+1):length(sig.list)) {
    overlap <- intersect(sig.list[[i]][[column_name]], sig.list[[j]][[column_name]])
    overlaps <- c(overlaps, overlap)
    percentage_overlap <- length(overlap) / length(unique(c(sig.list[[i]][[column_name]], sig.list[[j]][[column_name]]))) * 100
    percentage_overlaps <- c(percentage_overlaps, percentage_overlap)
  }
}


mean_overlap <- mean(percentage_overlaps)
cat(paste0("Mean percentage of overlap between ", column_name, " for all pairs of tables: ", round(mean_overlap, 2), "%\n"))

ovs <- table(overlaps)[order(table(overlaps), decreasing = T)][1:20]
ovs


ann <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/Annotation/gencode.v26.annotation.fixed.gtf.gz")
head(ann)

ann[ann$gene_id %in% names(ovs),]

sig.list.name <- c()

# Add gene name column:
for (i in 1:length(sig.list)) {
  sig.list.name[[i]] <- merge(sig.list[[i]], ann[,c("gene_id","gene_name")], by.x = "name", by.y = "gene_id")
}

sig.list.name[[1]]



all_genes <- unlist(lapply(sig.list.name, function(tbl) tbl$gene_name))
gene_counts <- table(all_genes)
gene_counts_df <- as.data.frame(gene_counts, stringsAsFactors = FALSE)
names(gene_counts_df) <- c("Gene_Name", "Occurrences")

gene_counts_df <- gene_counts_df %>%
  arrange(desc(Occurrences))


# We start selecting the genes that appear 5 or more times with ASE in our 12 samples:
gene_counts_df[gene_counts_df$Occurrences >= 5,]




twas <- readxl::read_xls("C://Users/Gerard/Downloads/41588_2023_1510_MOESM4_ESM.xls", sheet = 6, skip = 4)
head(twas)

twas.genes <- twas$Gene_name...2

intersect(twas.genes, gene_counts_df[gene_counts_df$Occurrences >= 6,]$Gene_Name)




# Load significant by diameter:
size <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Diameter.xlsx")
size$Gene.Name

column_name <- "gene_name"
percentage_overlaps <- c()
overlaps <- c()

for (i in 1:(length(sig.list.name))) {
  overlap <-
    intersect(sig.list.name[[i]][[column_name]], size$Gene.Name)
  overlaps <- c(overlaps, overlap)
  percentage_overlap <-
    length(overlap) / length(size$Gene.Name) * 100
  percentage_overlaps <-
    c(percentage_overlaps, percentage_overlap)
}

write.xlsx(table(overlaps)[order(table(overlaps), decreasing = T)], 
           "C://Users/Gerard/Desktop/AAA/RNAseq/ASE/Overlap_ASE_Diameter.xlsx", sep = "\t", colNames = F)

table(overlaps)[order(table(overlaps), decreasing = T)]
mean(percentage_overlaps)



