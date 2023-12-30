### Compare with previous results (Tunica-Specific):

library(openxlsx)

# Load our results:
cc <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Cases_Controls_No_Smoking.xlsx")
nrow(cc)

med <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 1)
med$Gene.Symbol <- gsub("-", "", med$Gene.Symbol)

#length(med$Gene.Symbol) # 5,889
#length(unique(med$Gene.Symbol)) # 5,507

#length(intersect(med$Gene.Symbol, cc$V1)) # 2,952

#length(setdiff(med$Gene.Symbol, cc$V1)) # 2,555
#length(setdiff(cc$V1, med$Gene.Symbol)) # 2,840

adv <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 2)
adv$Gene.Symbol <- gsub("-", "", adv$Gene.Symbol)


med.adv <- unique(c(med$Gene.Symbol, adv$Gene.Symbol))


# New genes compared to the previous larger study:
length(setdiff(cc$Gene.Name, med.adv)) # 3568.



#length(adv$Gene.Symbol) # 2,701
#length(unique(adv$Gene.Symbol)) #2,508

#length(intersect(adv$Gene.Symbol, cc$V1)) # 1,007

#length(setdiff(adv$Gene.Symbol, cc$V1)) # 1,501
#length(setdiff(cc$V1, adv$Gene.Symbol)) # 4,785


updown <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/UpDown Results.xlsx")

smla <- readxl::read_xls("C://Users/Gerard/Desktop/AAA/RNAseq/SmallLargeResults.xls")

rup <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/Ruptured.xlsx")

affyilu <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/Affyilu_Results.xlsx")

neck <- readxl::read_xls("C://Users/Gerard/Desktop/AAA/RNAseq/Neck.xls")

head(cc)
genes <- as.data.frame(cc[,c("Gene.Name","Beta")])


head(med , 1)
head(adv, 1)
head(updown, 1)
head(smla, 1)
head(rup, 1)
head(affyilu, 1)
head(neck, 1)


# Unifying gene lists
df_list <- list(adv = list(df = adv, code = "32907367"),
                med = list(df = med, code = "32907367"),
                updown = list(df = updown, code = "19111481"),
                smla = list(df = smla, code = "25944698"),
                rup = list(df = rup, code = "29191809"),
                affyilu = list(df = affyilu, code = "17634102"),
                neck = list(df = neck, code = "24529146"))

# Create new dataframe:
gene_presence <- data.frame(Gene.Name = genes$Gene.Name, Beta = genes$Beta, Symbol_Code = as.character(rep("", length(genes$Gene.Name))))
head(gene_presence)

# Check gene presence:
for (df_name in names(df_list)) {
  df_info <- df_list[[df_name]]
  in_gene_index <- gene_presence$Gene.Name %in% df_info$df$Gene.Symbol
  gene_presence[in_gene_index, "Symbol_Code"] <- ifelse(gene_presence[in_gene_index, "Symbol_Code"] == "", 
                                                        df_info$code,
                                                        paste(gene_presence[in_gene_index, "Symbol_Code"], df_info$code, sep = ","))
}

names(gene_presence) <- c("Gene Name","Beta","Previously Described")
gene_presence$'Functional Studies' <- NA


write.xlsx(gene_presence, "C://Users/Gerard/Desktop/AAA/RNAseq/AllGenes.xlsx")


new.genes <- gene_presence[gene_presence$`Previously Described` == "",]


write.xlsx(new.genes, "C://Users/Gerard/Desktop/AAA/RNAseq/NewGenes.xlsx")




