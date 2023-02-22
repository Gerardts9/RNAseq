library(data.table);library(openxlsx)

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



# In average, there are 530 genes with allelic specific expression in the 12 individuals with genotype information.



# Study overlap among significant genes between individuals:

sig.list <- lapply(groups, subset, binom_q < 0.05)

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

ann[ann$gene_id %in% names(ovs),]


ann <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/Annotation/gencode.v26.annotation.fixed.gtf.gz")
head(ann)
sig.list.name <- c()

# Add gene name column:
for (i in 1:length(sig.list)) {
  sig.list.name[[i]] <- merge(sig.list[[i]], ann[,c("gene_id","gene_name")], by.x = "name", by.y = "gene_id")
}

# Load significant by size:
size <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/GitHub/TopTable_Diameter.xlsx")
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

table(overlaps)[order(table(overlaps), decreasing = T)]
mean(percentage_overlaps)


# Load significant by type of aneurysm:
type <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/GitHub/TopTable_Type.xlsx")
column_name <- "gene_name"
percentage_overlaps <- c()

for (i in 1:(length(sig.list.name))) {
  overlap <-
    intersect(sig.list.name[[i]][[column_name]], type$Gene.Name)
  percentage_overlap <-
    length(overlap) / length(size$Gene.Name) * 100
  percentage_overlaps <-
    c(percentage_overlaps, percentage_overlap)
}

percentage_overlaps
mean(percentage_overlaps)

