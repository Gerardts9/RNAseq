library(dplyr);library(tidyr);library(openxlsx)

gwas.res <- tibble(readxl::read_xls("C://Users/Gerard/Desktop/AAA/RNAseq/AAAgen Supplementary.xls", sheet = 14, skip = 1))
class(gwas.res)
names(gwas.res)[names(gwas.res) == "Prioritized gene"] <- "Prioritized_Gene"


# Compare DEGs with GWAS:
all <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/AllGenes.xlsx")
class(all)
class(gwas.res)

merge_genes <- function(df, start_col, end_col) {
  df %>%
    rowwise() %>%
    mutate(Candidate_Genes = paste(
      unique(
        na.omit(
          unlist(
            strsplit(
              paste(na.omit(c_across(start_col:end_col)), collapse = ","),
              ","
            )
          )
        )
      ),
      collapse = ","
    )) %>%
    ungroup() %>%
    dplyr::select(1:3, "Prioritized_Gene", Candidate_Genes)
}


result <- as.data.frame(merge_genes(gwas.res, 11, 28))
head(result)


filter_genes <- function(df, gene_list) {
  df %>%
    rowwise() %>%
    mutate(DEGs_AAA_Control = paste(
      intersect(
        unlist(strsplit(Candidate_Genes, ",")),
        gene_list
      ),
      collapse = ","
    )) %>%
    ungroup()
}

result_with_matches <- filter_genes(result, all$Gene.Name)
head(result_with_matches)

nrow(result_with_matches[result_with_matches$DEGs_AAA_Control != "",])


head(result_with_matches)

write.xlsx(result_with_matches, "C://Users/Gerard/Desktop/AAA/RNAseq/DEGs_Candidates_Overlap.xlsx")


library(stringr)


result_filtered <- result_with_matches %>%
  filter(
    DEGs_AAA_Control != "" &                            # DEGs_AAA_Control is not empty
      !str_detect(DEGs_AAA_Control, Prioritized_Gene) &    # Prioritized gene is not in DEGs_AAA_Control
      str_detect(Candidate_Genes, Prioritized_Gene)        # Prioritized gene is in Candidate_Genes
  )

# Print the filtered result
head(result_filtered)

write.xlsx(result_filtered, "C://Users/Gerard/Desktop/AAA/RNAseq/DEGs_Candidates_Overlap_Changes.xlsx")




