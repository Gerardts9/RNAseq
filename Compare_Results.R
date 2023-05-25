library(openxlsx);library(readxl);library(data.table);library(clusterProfiler);library(org.Hs.eg.db);library(ggplot2)

ann <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/Annotation/gencode.v26.annotation.fixed.gtf.gz")
ann$gene_id <- substr(ann$gene_id, 1, 15)
head(ann)

# New results:
cc <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Cases_Controls.txt")

alt <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/SUPPA/AAA_Sig_Results.txt")
alt$gene_id <- substr(rownames(alt), 1, 15)
alt <- merge(alt, ann, by = "gene_id")


# Microarray studies:
med <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 1)
med$Gene.Symbol <- gsub("-", "", med$Gene.Symbol)

adv <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 2)
adv$Gene.Symbol <- gsub("-", "", adv$Gene.Symbol)

up <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/UpDown Results.xlsx", sheet = 1)
down <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/UpDown Results.xlsx", sheet = 2)

small <- read_excel("C://Users/Gerard/Desktop/AAA/RNAseq/SmallResults.xls")
large <- read_excel("C://Users/Gerard/Desktop/AAA/RNAseq/LargeResults.xls")

cc <- cc[!cc$V1 %in% med$Gene.Symbol,]
cc <- cc[!cc$V1 %in% adv$Gene.Symbol,]

cc <- cc[!cc$V1 %in% up$X2,]
cc <- cc[!cc$V1 %in% down$Abbreviation,]

cc <- cc[!cc$V1 %in% small$Symbol,]
cc <- cc[!cc$V1 %in% large$Symbol,]

#write.xlsx(cc[1:20,], "C://Users/Gerard/Desktop/ESHG/Poster/Table1V2.xlsx", sep = "\t")

common <- intersect(cc$V1, alt$gene_name)

#write.xlsx(cc[cc$V1 %in% common, ][1:10], "C://Users/Gerard/Desktop/ESHG/Poster/Table2.xlsx", sep = "\t")

gene_list <- unique(substr(cc$gene_id, 1, 15))

entrez_ids <- bitr(gene_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO enrichment analysis
go_enrich <- enrichGO(gene          = entrez_ids$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL", # Choose from "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)

# GO enrichment dot plot:
#png("C://Users/Gerard/Desktop/GO_Enrichment_New_Results.png", width = 4*800, height = 4*600, res = 300)
dotplot(go_enrich, showCategory = 10) + ggtitle("GO Enrichment")
#dev.off()

gene_id_list <- lapply(go_enrich@result$geneID, function(x) strsplit(x, "/")[[1]])

gene_symbols_df <- bitr(unlist(gene_id_list), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

gene_symbol_list <- lapply(gene_id_list, function(x) gene_symbols_df$SYMBOL[match(x, gene_symbols_df$ENTREZID)])
names(gene_symbol_list) <- names(gene_id_list)

collapsed_gene_symbols <- sapply(gene_symbol_list, function(x) paste(x, collapse = ","))
collapsed_gene_symbols_df <- data.frame(GeneSymbols = collapsed_gene_symbols, stringsAsFactors = FALSE)

df <- cbind(go_enrich@result, collapsed_gene_symbols_df)

df
