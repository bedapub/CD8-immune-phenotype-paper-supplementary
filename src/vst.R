train_vst <- function (matr, nsub = 1000, fitType = "parametric")
{
  if (nrow(matr) < nsub) {
    stop("less than 'nsub' rows,\n  it is recommended to use varianceStabilizingTransformation directly")
  }
  if (is.null(colnames(matr))) {
    colnames(matr) <- seq_len(ncol(matr))
  }
  object <- DESeq2::DESeqDataSetFromMatrix(
    matr,
    data.frame(row.names = colnames(matr)),
    ~1
  )
  object <- DESeq2::estimateSizeFactors(object)
  baseMean <- rowMeans(counts(object, normalized = TRUE))
  if (sum(baseMean > 5) < nsub) {
    stop("less than 'nsub' rows with mean normalized count > 5, \n  it is recommended to use varianceStabilizingTransformation directly")
  }
  object.sub <- object[baseMean > 5, ]
  baseMean <- baseMean[baseMean > 5]
  o <- order(baseMean)
  idx <- o[round(seq(from = 1, to = length(o), length = nsub))]
  object.sub <- object.sub[idx, ]
  object.sub <- DESeq2::estimateDispersionsGeneEst(object.sub, quiet = TRUE)
  object.sub <- DESeq2::estimateDispersionsFit(object.sub, fitType = fitType,
    quiet = TRUE)
  return(object.sub)

}

apply_vst <- function(matr, vst) {
  if (is.null(colnames(matr))) {
    colnames(matr) <- seq_len(ncol(matr))
  }
  object <- DESeq2::DESeqDataSetFromMatrix(
    matr,
    data.frame(row.names = colnames(matr)),
    ~1
  )
  object <- DESeq2::estimateSizeFactors(object)
  old_sf <- DESeq2::sizeFactors(vst)
  new_sf <- DESeq2::sizeFactors(object)
  samples_in_old_sf <- intersect(names(old_sf), names(new_sf))
  new_sf[samples_in_old_sf] <- old_sf[samples_in_old_sf]
  DESeq2::sizeFactors(object) <- new_sf
  suppressMessages({
    DESeq2::dispersionFunction(object) <- DESeq2::dispersionFunction(vst)
  })

  vsd <- DESeq2:: varianceStabilizingTransformation(object, blind = FALSE)
  return(SummarizedExperiment::assay(vsd))
}
