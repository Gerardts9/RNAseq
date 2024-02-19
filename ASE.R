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
library(dplyr);library(biomaRt);library(aPEAR);library(factoextra)
library(clusterProfiler);library(org.Hs.eg.db)


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
cases <- as.data.frame(fread("/home/gerard/AAA/ASE/Cases_ratio.txt"))
cases <- cases %>% remove_rownames %>% column_to_rownames(var="V1")

controls <- as.data.frame(fread("/home/gerard/AAA/ASE/Controls_ratio.txt"))
controls <- controls %>% remove_rownames %>% column_to_rownames(var="V1")

common <- intersect(names(cases),names(controls))

cases <- cases[common]
controls <- controls[common]

if (identical(names(cases),names(controls)) != TRUE){warning("Genes are not the same")}

all <- rbind(cases,controls)
all[1:5,1:5]
dim(all) # 399 56200

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



# Compare p-values between candidate genes and all genes:
######################
# Load AAAGen results:

gwas <- fread("/home/gerard/AAA_MR/AAA_hg38.txt.gz")
gwas2 <- gwas[gwas$'P-value' < 5e-8, ]


# Concatenate all files:
aes <- list.files("/home/gerard/AAA/ASE/Results/", pattern = "allelic_counts.txt", full.names = T, recursive = T)
over.8 <- c()

for (i in 1:length(aes)) {
  ff <- fread(aes[i])
  ff$Table <- i
  ff2 <- subset(ff, ff$totalCount >= 8)
  over.8 <- rbind(over.8, ff2)
}

# Perform a binomial test for deviation from 0.5:
over.8$binom_p <- apply(over.8[,c("refCount","altCount")], 1, function(x) binom.test(x[1],x[1]+x[2],p=0.5)$p.value)
over.8$contig <- gsub("chr","",over.8$contig)
over.8$MarkerName <- paste0(over.8$contig, ":", over.8$position)

length(intersect(gwas2$MarkerName, over.8$MarkerName))

over.8[over.8$binom_p == 0,]$binom_p <- 5e-324

png("/home/gerard/AAA/ASE/Plots/QQPlot_All_Genes.png")
ggd.qqplot(over.8$binom_p)
dev.off()

png("/home/gerard/AAA/ASE/Plots/QQPlot_Candidate_Genes.png")
ggd.qqplot(over.8[over.8$MarkerName %in% gwas2$MarkerName,]$binom_p)
dev.off()


png("/home/gerard/AAA/ASE/Plots/QQPlot.png")
o1 = -log10(sort(over.8$binom_p, decreasing = FALSE))
e1 = -log10(1:length(o1) / length(o1))

o2 = -log10(sort(over.8[over.8$MarkerName %in% gwas2$MarkerName,]$binom_p, decreasing = FALSE))
e2 = -log10(1:length(o2) / length(o2))

max_e = max(max(e1), max(e2))
max_o = max(max(o1), max(o2))

plot(e1, o1, pch = 19, cex = 1, 
     xlab = expression(Expected ~~-log[10](italic(p))), 
     ylab = expression(Observed ~~-log[10](italic(p))), 
     xlim = c(0, max_e), ylim = c(0, max_o), col = "red")

points(e2, o2, pch = 19, cex = 1, col = "blue")
dev.off()
######################



# Prepare ratios counts on cases and controls:
######################
cases <- fread("/home/gerard/AAA/ASE/Cases_Counts.txt")
cases <- cases %>% remove_rownames %>% column_to_rownames(var="V1")

controls <- fread("/home/gerard/AAA/ASE/Controls_Counts.txt")
controls <- controls %>% remove_rownames %>% column_to_rownames(var="V1")

common <- intersect(names(cases),names(controls))

cases <- cases[common]
controls <- controls[common]

if (identical(names(cases),names(controls)) != TRUE){warning("Genes are not the same")}


# Find columns containing "_aCount":
aCount_cols <- grep("_aCount", colnames(cases))

# Extract gene names from column names:
gene_names <- sub("_aCount", "", colnames(cases)[aCount_cols])

# Loop through each gene and calculate the total counts for CASES:

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


# Calculate ratio between aCount and totalCount for CASES:

for (gene in gene_names) {
  # Select columns for aCount and totalCount for the current gene:
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

#fwrite(cases, "/home/gerard/AAA/ASE/Cases_ratio.txt", sep = "\t", row.names = TRUE)




# Loop through each gene and calculate the sum of counts for CONTROLS:

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


# Calculate ratio between aCount and totalCount for CONTROLS:

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

#fwrite(controls, "/home/gerard/AAA/ASE/Controls_ratio.txt", sep = "\t", row.names = TRUE)
######################



# Compare ratios between cases and controls and save results from Wilcoxon non-parametric test:
######################
cases <- fread("/home/gerard/AAA/ASE/Cases_ratio.txt")
cases <- cases %>% remove_rownames %>% column_to_rownames(var="V1")

controls <- fread("/home/gerard/AAA/ASE/Controls_ratio.txt")
controls <- controls %>% remove_rownames %>% column_to_rownames(var="V1")

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
p_values[is.na(p_values)] <- 1 # Set p-values equal to NA to 1.


p_values <- data.frame(p_values)
p_values$Gene_ID <- rownames(p_values)
names(p_values) <- c("P.value","Gene")
p_values <- p_values[,c("Gene","P.value")]

p_values$FDR <- p.adjust(p_values$P.value, method = "fdr")


# Save results:
#write.table(p_values, "/home/gerard/AAA/ASE/Pvalues.txt", sep = "\t", col.names = T, quote = F, row.names = F)
######################



# Identify significant ASE patterns and perform cluster analysis:
######################
p_values <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/ASE/Pvalues.txt")
head(p_values)

sig.genes <- p_values[p_values$FDR < 0.05]

nrow(sig.genes) # 1815


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_ids <- substr(sig.genes$Gene, 1, 15)

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

set.seed(28)
tiff("C://Users/Gerard/Desktop/AAA/RNAseq/ASE/ASE_Cluster4.tiff", height = 3500, width = 3500, res = 250)
plotPathClusters(
  enrichment = go_enrich@result,
  sim = clusters$similarity,
  clusters = clusters$clusters,
  fontSize = 10,
  outerCutoff = 0.1,
  drawEllipses = TRUE,
  nodeSize = "Count",
  colorBy = "pvalue",
  colorType = 'pval',
  repelLabels = TRUE
)
dev.off()
######################



# Load allelic counts:
######################
alelic.counts.files <- list.files("/home/gerard/AAA/ASE/Results/", pattern = "allelic_counts.txt", full.names = T, recursive = T)
alelic.counts <- c()

for (i in 1:length(alelic.counts.files)) {
  ff <- fread(alelic.counts.files[i])
  ff$Table <- i
  ff2 <- subset(ff, ff$totalCount >= 8)
  alelic.counts <- rbind(alelic.counts, ff2)
}

# Perform a binomial test for deviation from 0.5:
alelic.counts$binom_p <- apply(alelic.counts[,c("refCount","altCount")], 1, function(x) binom.test(x[1],x[1]+x[2],p=0.5)$p.value)

# Perform multiple testing correction with FDR for the 12 individuals:
alelic.counts$binom_q <- p.adjust(alelic.counts$binom_p, method = "fdr")
######################



# Load gene counts:
######################
# Concatenate all files:
gene.counts.files <- list.files("/home/gerard/AAA/ASE/Results/", pattern = "case_gene_ae.txt", full.names = T, recursive = T)
gene.counts <- c()

for (i in 1:length(gene.counts.files)) {
  ff <- fread(gene.counts.files[i])
  ff$Table <- i
  ff2 <- subset(ff, ff$totalCount >= 8)
  gene.counts <- rbind(gene.counts, ff2)
}

# Perform a binomial test for deviation from 0.5:
gene.counts$binom_p <- apply(gene.counts[,c("aCount","bCount")], 1, function(x) binom.test(x[1],x[1]+x[2],p=0.5)$p.value)

# Perform multiple testing correction with FDR for the 12 individuals:
gene.counts$binom_q <- p.adjust(gene.counts$binom_p, method = "fdr")
######################


# Load GWAS summary statistics:
######################
# Prepare GWAS summary statistics:
gwas <- fread("/home/gerard/AAA_MR/AAA_hg38.txt.gz")
#gwas <- gwas[gwas$'P-value' < 5e-8, ]
#fwrite(gwas, "/home/gerard/AAA_MR/AAA_hg38_sig.txt.gz", sep = "\t")
#gwas <- fread("/home/gerard/AAA_MR/AAA_hg38_sig.txt.gz")
######################




# Obtain number of ASE genes:

groups <- split(gene.counts, gene.counts$Table)

sig.rows <- function(df, columna) {
  sum(df[[columna]] < 0.05, na.rm = TRUE)
}

sig.rows.table <- sapply(groups, function(df) {
  sig.rows(df, "binom_q")
})

mean(sig.rows.table) # 529.17

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


ann <- fread("/home/gerard/AAA/refs/gencode.v26.annotation.fixed.gtf.gz")
head(ann)

sig.list.name <- c()

# Add gene name column:
for (i in 1:length(sig.list)) {
  sig.list.name[[i]] <- merge(sig.list[[i]], ann[,c("gene_id","gene_name")], by.x = "name", by.y = "gene_id")
}

sig.list.name[[1]]

gene_names <- unlist(lapply(sig.list.name, function(tbl) tbl$gene_name))
ensembl_names <- unlist(lapply(sig.list.name, function(tbl) tbl$name))

all_genes <- data.frame(Gene_Name = gene_names, Gene_ID = ensembl_names)
head(all_genes)

gene_counts <- table(all_genes$Gene_Name)
gene_counts_df <- as.data.frame(gene_counts, stringsAsFactors = FALSE)
names(gene_counts_df) <- c("Gene_Name", "Occurrences")

all_genes <- unique(all_genes)
gene_counts_df <- merge(gene_counts_df, all_genes, by = "Gene_Name")

gene_counts_df <- gene_counts_df %>%
  arrange(desc(Occurrences))

head(gene_counts_df)


# We select the genes that appear 5 or more times with ASE in our 12 samples:

sel.genes <- gene_counts_df[gene_counts_df$Occurrences >= 5,]
nrow(sel.genes) # 90
sel.genes


# Load eQTLs in blood and artery aorta:
eqtl.blood <- fread("/home/gerard/Colocalization_Gtex/Signif/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz")
eqtl.aorta <- fread("/home/gerard/Colocalization_Gtex/Signif/Artery_Aorta.v8.signif_variant_gene_pairs.txt.gz")

eqtl.blood$variant_id <- gsub("_b38", "", eqtl.blood$variant_id)
eqtl.aorta$variant_id <- gsub("_b38", "", eqtl.aorta$variant_id)

# Compare with genes that have different patterns between AAA cases and GTEx controls:
p_values <- fread("/home/gerard/AAA/ASE/Pvalues.txt")

sig.ase.genes <- p_values[p_values$FDR < 0.05]

nrow(sig.ase.genes) # 1815
head(sig.ase.genes)

gene_counts_df[gene_counts_df$Gene_ID %in% intersect(sel.genes$Gene_ID, sig.ase.genes$Gene),] 
# SNURF gene



# Compare if genes are also present in the TWAS of AAAgen:
twas <- readxl::read_xls("/home/gerard/AAA/ASE/AAAgen Supplementary.xls", sheet = 6, skip = 4)
head(twas)

twas.genes <- twas$Gene_name...2

intersect(sel.genes$Gene_Name, twas.genes)
# THBS2 gene.


# Compare if genes are also present in the GWAS of AAAgen:
gwas.res <- readxl::read_xls("/home/gerard/AAA/ASE/AAAgen Supplementary.xls", sheet = 14, skip = 1)
head(gwas.res)

intersect(sel.genes$Gene_Name, gwas.res$`Prioritized gene`)
# SPP1 and THBS2 genes.


#########
# SNURF #
#########


# Extract genetic variants in SNURF genes:
sig.df <- do.call(rbind, sig.list.name)
head(sig.df)
snurf.vars <- unique(unlist(strsplit(as.character(sig.df[sig.df$gene_name %in% "SNURF",]$variants), ",")))
snurf.vars

# Check variants in allelic counts: We see that chr15_24974365_T_C SNP determines expression in 8 samples.
alelic.counts[alelic.counts$variantID %in% snurf.vars &
                alelic.counts$binom_q < 0.05,]


sum(alelic.counts[alelic.counts$variantID == "chr15_24974365_T_C",]$refCount)
sum(alelic.counts[alelic.counts$variantID == "chr15_24974365_T_C",]$altCount)
# We see that T allele of chr15_24974365_T_C SNP is overexpressed than C allele.


# Check variants in GWAS: Not significant locus.
min(gwas[gwas$chr == 15 & gwas$pos > 24966466 - 500000 & gwas$pos < 24974365 + 500000,]$'P-value')
gwas[gwas$chr == 15 & gwas$pos == 24974365,]


# Chek eQTL:
eqtl.blood[grep("chr15_24974365", eqtl.blood$variant_id),] # In Blood, eQTL of SNRPN.
eqtl.aorta[grep("chr15_24974365", eqtl.aorta$variant_id),] # In Artery Aorta, eQTL of lnc-SNRPN-8.



#########
# THBS2 #
#########

thbs.vars <- unique(unlist(strsplit(as.character(sig.df[sig.df$gene_name %in% "THBS2",]$variants), ",")))
thbs.vars

thbs.counts <- alelic.counts[alelic.counts$variantID %in% thbs.vars &
                               alelic.counts$binom_q < 0.05, ]
# There are many variants that determine the expression of THBS2 gene.

thbs.vars.pos <- str_split_fixed(thbs.vars, "_",4)[,2]

# Check gwas:
thbs.gwas <- gwas[gwas$chr == 6 & gwas$pos %in% thbs.vars.pos & gwas$'P-value' < 5e-8,]

common.positions <- intersect(thbs.counts$position, thbs.gwas$pos)

thbs.counts[thbs.counts$position %in% common.positions,]
thbs.gwas[thbs.gwas$pos %in% common.positions,]



# Check eQTL of TopSNP: 6:169177891: # No significant eQTL.

# Check for common SNPs between GWAS and ASE:
eqtl.blood[grep("chr6_169220665", eqtl.blood$variant_id),] # No significant eQTL.
eqtl.aorta[grep("chr6_169220665", eqtl.aorta$variant_id),] # eQTL of THBS2 gene.

eqtl.blood[grep("chr6_169222395", eqtl.blood$variant_id),] # No significant eQTL.
eqtl.aorta[grep("chr6_169222395", eqtl.aorta$variant_id),] # eQTL of THBS2 gene.


########
# SPP1 #
########


# Extract genetic variants in SNURF genes:
sig.df <- do.call(rbind, sig.list.name)
head(sig.df)
spp.vars <- unique(unlist(strsplit(as.character(sig.df[sig.df$gene_name %in% "SPP1",]$variants), ",")))
spp.vars

spp.counts <- alelic.counts[alelic.counts$variantID %in% spp.vars &
                               alelic.counts$binom_q < 0.05, ]
# There are 6 variants that determine the expression of SPP1 gene.
table(spp.counts$variant)


spp.vars.pos <- str_split_fixed(spp.vars, "_",4)[,2]

# Check gwas:
spp.gwas <- gwas[gwas$chr == 4 & gwas$pos %in% spp.vars.pos & gwas$'P-value' < 5e-8,]

common.positions <- intersect(spp.counts$position, spp.gwas$pos)

spp.counts[spp.counts$position %in% common.positions,]
spp.gwas[spp.gwas$pos %in% common.positions,]



# Check eQTL of TopSNP:
eqtl.blood[grep("chr4_87854008", eqtl.blood$variant_id),] # eQTL of SPP1 and PKD2 genes.
eqtl.aorta[grep("chr4_87854008", eqtl.aorta$variant_id),] # eQTL of SPP1 gene.



# Check for common SNPs between GWAS and ASE:
eqtl.blood[grep("chr6_169220665", eqtl.blood$variant_id),] #
eqtl.aorta[grep("chr6_169220665", eqtl.aorta$variant_id),] #







# Compare if genes are also DEGs between cases and controls, and by diameter:
cc <- read.xlsx("/home/gerard/AAA/SigRes/Cases_Controls_No_Smoking.xlsx")
dia <- read.xlsx("/home/gerard/AAA/SigRes/Diameter.xlsx")


intersect(sel.genes$Gene_Name, cc$Gene.Name) # 12 genes between selected and case-control.
intersect(sel.genes$Gene_Name, dia$Gene.Name) # No genes between selected and diameter.










