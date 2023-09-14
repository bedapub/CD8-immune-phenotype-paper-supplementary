#' Create a matrix with symbols used for row.names.
#' @param eset ExpressionSet
#' @return A matrix.
eset_symbol_matrix <- function(eset) {
  m <- eset %>%
      Biobase::exprs()
  rownames(m) <- if_else(
      is.na(Biobase::fData(eset)$symbol) | Biobase::fData(eset)$symbol == "",
      rownames(m),
      Biobase::fData(eset)$symbol
  ) %>%
    make.names(unique = TRUE)
  m
}

#' Generate data.frame for ML model
#' @param eset ExpressionSet
#' @param use_genes Include gene expression values.
#' @param exclude_genes Remove some genes, character vector.
#' Only makes sense if `use_genes` is `TRUE`.
#' @param include_only_genes Only keep these genes, character vector.
#' Only makes sense if `use_genes` is `TRUE`.
#' @param signatures A matrix with signature scores.
#' @param signature_namespace Include signatures starting form prefix + `"_"` symbol.
create_modeling_df <- function(eset,
                               use_genes = TRUE,
                               exclude_genes = NULL,
                               include_only_genes = NULL,
                               signatures = NULL,
                               signature_namespace = NULL) {
  if ("CD8IMMPH" %in% colnames(Biobase::pData(eset))) {
    df <- data.frame(IHC_Status = str_to_lower(eset$CD8IMMPH))
  } else {
    df <- data.frame(IHC_Status = rep(NA_character_, ncol(eset)))
  }
  rownames(df) <- colnames(eset)
  df$IHC_Status <- factor(df$IHC_Status,
                          levels = c("desert", "excluded", "inflamed"))
  if (use_genes) {
    m_genes <- eset_symbol_matrix(eset) %>%
      t()
    if (!is.null(exclude_genes)) {
      m_genes <- m_genes[, !colnames(m_genes) %in% exclude_genes]
    }
    if (!is.null(include_only_genes)) {
      m_genes <- m_genes[, colnames(m_genes) %in% include_only_genes]
    }
    df <- cbind(df, m_genes)
  }
  if (!is.null(signatures)) {
    if (is.null(signature_namespace)) {
      warning("signature_namespace is NULL, using all the signatures")
    } else {
      signatures <- signatures[str_starts(rownames(signatures), str_c(signature_namespace, "_")), ]
    }
  }
  if (!is.null(signatures)) {
    signatures <- t(signatures[, rownames(df)])
    stopifnot(nrow(signatures) == nrow(df))
    stopifnot(all(rownames(signatures) == rownames(df)))
    df <- cbind(df, signatures)
  }
  df
}


#' Matches cit_ext harmonized location to Atezo
#' @param eset citext_anno_BMFLOC the BMLOC harmonized from cit_ext
#' @param citmain_anno_BMFLOC the BMLOC semi-harmonized from atezo
#' @return A vector to replace the BMLOC in Atezo
match_citext_atezo_loc <- function(citext_anno_BMFLOC,
                                   citmain_anno_BMFLOC) {

  cats= gsub("_", " ", unique(citext_anno_BMFLOC))

  tmp <- data.frame(BMFLOC=citmain_anno_BMFLOC)
  tmp$newBMFLOC <- NA

  # First easy fix
  for( i in 1:length(cats)) {
    x=cats[i]
    tmp$newBMFLOC [ grep (x, tmp$BMFLOC, ignore.case = TRUE) ] <- x
  }

  # Exceptions
  tmp <- tmp %>% dplyr::mutate(newBMFLOC = ifelse( grepl("ABDOM", .data$BMFLOC), "abdominal_cavity", .data$newBMFLOC )) %>%
                dplyr::mutate(newBMFLOC = ifelse( grepl("ABD WALL", .data$BMFLOC), "abdominal_cavity", .data$newBMFLOC ))  %>%
                dplyr::mutate(newBMFLOC = ifelse( grepl("BRONCHIC", .data$BMFLOC), "bronchus", .data$newBMFLOC )) %>%
                dplyr::mutate(newBMFLOC = ifelse( grepl("URET", .data$BMFLOC), "urothelial", .data$newBMFLOC )) %>%  # tbc if urether could be classified as urothelial
                dplyr::mutate(newBMFLOC = ifelse( grepl("LYMPH|LYMF|LYMPHNODE", .data$BMFLOC), "lymph_node", .data$newBMFLOC )) %>%
                dplyr::mutate(newBMFLOC = ifelse( grepl("PELVIC", .data$BMFLOC), "pelvis", .data$newBMFLOC )) %>%
                dplyr::mutate(newBMFLOC = ifelse( grepl("PLEURA", .data$BMFLOC), "pleura_and_mediastinum", .data$newBMFLOC ))

  # For those that fit nothing
  tmp$newBMFLOC [is.na(tmp$newBMFLOC )] = "other"

  # to revert the original cats "_"

  tmp$newBMFLOC = gsub(" ", "_", tmp$newBMFLOC)

  return(tmp$newBMFLOC)
}

#' Assigns whether primary or metastatic based on indication vs location.
#' It provides primary, as logical as possible (see few exceptions in function), and multiple metastatic states or as "other".
#' @param data_df dataframe that should contain BMFLOC and INDICAT as columns
#' @return A vector to add as STATUS

 assign_status <- function (data_df) {

  data_df$STATUS = NA

  # Defining primary vs Mets
  data_df$STATUS[  which( data_df$INDICAT == toupper(data_df$BMFLOC)) ] = "Primary"
  data_df$STATUS[  - which( data_df$INDICAT == toupper(data_df$BMFLOC))] = "Other Mets"

  # defining exceptions or other primary
  data_df <- data_df %>%  dplyr::mutate(STATUS = ifelse( .data$INDICAT == "COLON" & .data$BMFLOC == "abdominal_cavity", "Primary", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$INDICAT == "CERVIX UTERI" & .data$BMFLOC == "uterus_includesCervix", "Primary", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$INDICAT == "NSCLC" & .data$BMFLOC == "lung", "Primary", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$INDICAT == "SOFT TISSUE" & .data$BMFLOC == "soft_tissue", "Primary", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$INDICAT == "OVARY" & .data$BMFLOC == "abdominal_cavity", "Primary", .data$STATUS )) %>%   # This seems to be Primary peritoneal cancer (PPC) the most common ovarian cancer cells
                          dplyr::mutate(STATUS = ifelse( .data$INDICAT == "STOMACH" & .data$BMFLOC == "abdominal_cavity", "Primary", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$INDICAT == "UROTHELIAL BLADDER CANCER" & .data$BMFLOC == "bladder", "Primary", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$INDICAT == "UROTHELIAL BLADDER CANCER" & .data$BMFLOC == "urothelial", "Primary", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$INDICAT == "PLEURA" & .data$BMFLOC == "pleura_and_mediastinum", "Primary", .data$STATUS ))

  # Adding some mets
  data_df <- data_df %>%  dplyr::mutate(STATUS = ifelse( .data$STATUS == "Other Mets" & .data$BMFLOC == "liver", "Liver Mets", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$STATUS == "Other Mets" & .data$BMFLOC == "lung", "Lung Mets", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$STATUS == "Other Mets" & .data$BMFLOC == "lymph_node", "LN Mets", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$STATUS == "Other Mets" & .data$BMFLOC == "kidney", "Kidney Mets", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$STATUS == "Other Mets" & .data$BMFLOC == "head_and_neck", "Head/Neck Mets", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$STATUS == "Other Mets" & .data$BMFLOC == "breast", "Breast Mets", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$STATUS == "Other Mets" & .data$BMFLOC == "breast", "Breast Mets", .data$STATUS )) %>%
                          dplyr::mutate(STATUS = ifelse( .data$STATUS == "Other Mets" & .data$BMFLOC == "pancreas", "Pancreas Mets", .data$STATUS ))

   return(data_df$STATUS)
  }
