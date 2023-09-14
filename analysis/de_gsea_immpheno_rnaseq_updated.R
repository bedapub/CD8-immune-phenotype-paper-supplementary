#' ---
#' title: "RNAseq HGX Immune phenotypes DE and GSEA analysis"
#' subtitle: "RT squad - Immune phenotypes"
#' author: "Nicolas Staedler"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'     includes:
#'        in_header: !expr system.file("extdata", "header.html", package = "ecfc")
#'        after_body: !expr system.file("extdata", "footer.html", package = "ecfc")
#'     css: !expr system.file("extdata", "style.css", package = "ecfc")
#'     mathjax: default
#'     number_sections: yes
#'     toc: yes
#'     toc_depth: 3
#'     toc_float:
#'       collapsed: true
#'       smooth_scroll: true
#'     code_folding: hide
#'  pdf_document:
#'    number_sections: yes
#'    toc: yes
#'    toc_depth: 3
#' params:
#'  interactfl: FALSE
#'  output_collection: "UUID"
#'  data_collection: "UUID"
#' ---
#'

#+ results = 'hide', message = FALSE
# global chunk options
library(knitr)
opts_chunk$set(echo = TRUE, eval = TRUE, warning = FALSE, message = FALSE, fig.align = 'center', cache = FALSE)

# load packages
library(BioQC)
library(edgeR)
library(aws.s3)
library(tidyverse)
library(ecfc)
library(ComplexHeatmap)
library(DT)
library(GSEABase)
library(ggrepel)
library(ggpubr)
library(plotly)
library(RColorBrewer)
library(openxlsx)
library(purrr)
library(UpSetR)
library(VennDiagram)

#' ***
#' # Content
#'
#' Identify genes and genesets which are differentially expressed
#'
#' - desert vs inflamed
#' - excluded vs inflamed
#' - excluded vs desert.
#'
#' We add TISSUE as a covariate.
#'

#' ***
#' # Results
#'
#+ include=FALSE
#params <- list(output_collection="UUID")

#+
# use Iakov's expression data
eset.counts <- aws.s3::s3readRDS(
  "counts_eset.Rds",
  bucket = params$output_collection
)

eset.counts <- eset.counts[, !is.na(eset.counts$CD8IMMPH)]

#+
counts <- exprs(eset.counts)
samples <- pData(eset.counts)
bmfloc_nsel <- names(table(samples$BMFLOC))[table(samples$BMFLOC)>=100]
samples$TISSUE <- ifelse(samples$BMFLOC%in%bmfloc_nsel,samples$BMFLOC,"OTHER_UNK")
genes <- fData(eset.counts)
group <- samples$CD8IMMPH
dge <- DGEList(counts=counts,
               samples=samples,
               genes = genes,
               group=group)
rownames(dge$counts) <- dge$genes$symbol

# filtering
keep.exprs <- filterByExpr(dge)
dge2 <- dge[keep.exprs,, keep.lib.sizes=FALSE]

# normalization using TMM (library size --> effective library size)
dge3 <- calcNormFactors(dge2, method = "TMM")

# make expression set (log cpm and log rpkm)
pdat <- dge3$samples
emat <- cpm(dge3,log=TRUE) #normalized log cpm
emat.rpkm <- rpkm(dge2, gene.length="size", log=TRUE, normalized.lib.sizes=TRUE)
annot <- dge3$genes
rownames(annot) <- annot$GeneSymbol
rownames(emat) <- rownames(annot)
rownames(emat.rpkm) <- rownames(annot)
eset <- ExpressionSet(assayData=emat)
pData(eset) <- pdat
fData(eset) <- annot
eset.rpkm <- ExpressionSet(assayData=emat.rpkm)
pData(eset.rpkm) <- pdat
fData(eset.rpkm) <- annot


#+
# get signature scores
gsets <- list(
  sc = "CellNames_scseqCMs6_sigs_for_Cd8.gmt"
 #cancer = "c4.cm.v7.4.symbols.gmt"
 #  cancer = "c6.all.v7.4.symbols.gmt"
 #  hallmark = "h.all.v7.4.symbols.gmt"
) %>%
  imap(~ aws.s3::s3read_using(BioQC::readGmt, namespace = .y,
                              object = .x, bucket = params$data_collection)) %>%
  rlang::flatten() %>%
  BioQC::GmtList()


gsets <- gsets$sc

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


signature_scores <- BioQC::wmwTest(eset_symbol_matrix(eset),  ## STOP where is this function? otherwise use the gmts somewhere ...
                             gsets,
                             valType="r"
)


#'
#' ## Data overview {.tabset .tabset-pills}
#'
#' ### Zero count genes
#'
#'  Genes with zero counts for all samples.
#'
table(rowSums(dge$counts==0)==ncol(dge))

#' ### Low expressed genes
#'
#' Low expressed genes.
#+
table(!keep.exprs)

#' ### Density plot
#'
#' Densitiy plot before and after removing low expressed genes.
#'
#+
ss <- sample(1:ncol(dge),500)
x <- dge[,ss]
lcpm <- cpm(x, log=TRUE)
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(x)
col <- factor(x$samples$CD8IMMPH)
levels(col) <- c("green","blue","red")
col <- as.character(col)
par(mfrow=c(1,2))
plot(density(lcpm[,1]), ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y)
}

x <- dge2[,ss]
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y)
}

#'
#' ### Samples
#'
#' Overview of samples
#+
#+ fig.width = 8,fig.height = 6, results="asis"
dd <- pdat%>%
  dplyr::select(TISSUE,INDICAT,CD8IMMPH)

dd%>%
  ggplot(.,aes(x=TISSUE,fill=CD8IMMPH))+
  geom_bar(position = position_dodge())+
  theme(text=element_text(size=15),
        legend.position = "down")

#'
#' ## MDS plots {.tabset .tabset-pills}
#'
#' ### By CD8IMMPH
#+
#+ fig.width = 6,fig.height = 6, results="asis"
x <- dge3[,ss]
lcpm <- cpm(x, log=TRUE)
grp <- factor(x$samples$CD8IMMPH)
col.grp <- grp
levels(col.grp) <-  c("#7570b3","#1b9e77","#d95f02")
mds <- plotMDS(lcpm,gene.selection = "common",plot=FALSE)
plotMDS(mds,pch=3,col=as.character(col.grp),cex=0.75)
legend("topright", levels(grp), text.col=levels(col.grp), bty="n")

#' ### By TISSUE
#+
#+ fig.width = 6,fig.height = 6, results="asis"
grp <- factor(x$samples$TISSUE)
col.grp <- grp
levels(col.grp) <-  colorRampPalette(brewer.pal(8, "Set1"))(nlevels(grp))
plotMDS(mds,pch=3,col=as.character(col.grp),cex=0.75)
legend("topright", levels(grp), text.col=levels(col.grp), bty="n")

#' ### List of Genesets
#+
t.d.l <- lapply(gsets,function(x){
  t.name <- x$name
  t.desc <- x$desc
  t.genes <- paste(x$genes[x$genes!=""],collapse=", ")
  out <- data.frame(name=t.name,description=t.desc,genes=t.genes)
  return(out)
})
t.d <- do.call(rbind,t.d.l)
DT::datatable(t.d,rownames=FALSE,caption="List of gene sets")

#'
#' ## Differential expression {.tabset .tabset-pills}
#'
#+ fig.width = 6,fig.height = 6, results="asis"
# logFC cutoff used to highlight significant hits in the plot
volcano.logFC.cutoff <- log2(1.5)
# number of significant genes for which labels are shown on both sides (up- and down-regulated)
volcano.numLabels <- 20
# fdr cutoffs
alpha.sig <- 0.05
# number of top genes to plot
ntop <- 9
# topflag
topfl <- quo(pval.adj<alpha.sig&(abs(logFC)>volcano.logFC.cutoff))

#+
# design matrix and contrast
myform <- ~0+CD8IMMPH+TISSUE

design.x <-  model.matrix(myform,data=pdat)
colnames(design.x) <- gsub("CD8IMMPH","",colnames(design.x))
linfct <- makeContrasts(INFLAMED_DESERT=INFLAMED-DESERT,
                        INFLAMED_EXCLUDED=INFLAMED-EXCLUDED,
                        EXCLUDED_DESERT=EXCLUDED-DESERT,levels=design.x)

#+
# limma voom analysis
v <- voom(counts=dge3,design=design.x,plot=TRUE)
lmfit <- lmFit(v,design=design.x)
cfit <- contrasts.fit(lmfit,contrasts=linfct)
cfitbayes <- eBayes(cfit)

#+
# put results together
t.d <- data.frame(sigma=cfitbayes$sigma,
                  amean=cfitbayes$Amean,
                  PARAM=rownames(cfitbayes$coefficients))
res.l <- sapply(colnames(linfct),FUN=function(i){
  tt <- topTable(cfitbayes,coef=i,adjust.method = "BH",number=Inf)
  tt$CONTR <- i
  tt$PARAM <- tt$symbol
  return(tt)}, simplify = FALSE,USE.NAMES = TRUE
)
res <- data.frame(do.call("rbind",res.l))%>%
  dplyr::select(PARAM,CONTR,logFC,P.Value,adj.P.Val)%>%
  dplyr::rename(pval=P.Value,pval.adj=adj.P.Val)%>%
  dplyr::mutate(UPDOWN=ifelse(logFC>0,"UP","DOWN"))%>%
  dplyr::mutate(TOPFLAG=!!topfl)%>%
  group_by(CONTR,UPDOWN)%>%
  arrange(pval)%>%
  dplyr::mutate(TOPFLAG2=ifelse(between(row_number(), 1, volcano.numLabels),TRUE,FALSE))%>%
  dplyr::mutate(TOPFLAG2=TOPFLAG2&TOPFLAG)%>%
  ungroup
res1 <- left_join(res,t.d,by="PARAM")
res2 <- left_join(res1,
                  fData(eset),
                  by=c("PARAM"="symbol")
)

#' ### Mean-variance relationship
#+
plotSA(cfitbayes)

#+
# prepare volcano plot
mycontr <- unique(res2$CONTR)

ggpl <- lapply(mycontr,FUN=function(x){

  t.ggpl <- res2%>%
    dplyr::filter(CONTR==x)%>%
    dplyr::mutate(nLogPval=ifelse(pval==0,round(-log10(min(pval[pval>0]))),-log10(pval)))%>%
    dplyr::mutate(adjPval=round(pval.adj,4))%>%
    volcano_plot(.,
                 formula=nLogPval~logFC+UPDOWN,
                 group="PARAM",transform=NULL,topflag="TOPFLAG2")+
    labs(title = x)+
    theme(text=element_text(size=20))+
    xlim(c(-4,4))+ylim(c(0,150))

  return(t.ggpl)

}
)
names(ggpl) <- mycontr


#' ### Volcano plots {.tabset}
#+ results="asis"
#+ fig.width = 6,fig.height = 6, results="asis"
if(params$interactfl){
  htmltools::tagList(lapply(ggpl,FUN=
                              function(x){
                                x1 <- x%+%aes(label=logFC,label2=adjPval,label3=PARAM)+
                                  theme(legend.position = "top")

                                ggplotly(x1,tooltip=c("logFC","adjPval","PARAM"),
                                         height=700,width=800)
                              })
  )
}else{
  lapply(ggpl,
         function(x){
           t.ggpl <- x+
             geom_text_repel(aes(label=ifelse(TOPFLAG2,PARAM,"")),
                             size=3.5)
           t.ggpl

         })
}


ggpl <-lapply(ggpl,
       function(x){
         t.ggpl <- x+
           geom_text_repel(aes(label=ifelse(TOPFLAG2,PARAM,"")),
                           size=3.5)
       })
for(ct in names(ggpl)){
wd=getwd()
svg(paste(wd, "/volcano_",ct, ".svg", sep=""), width = 8, height = 8)
print(ggpl[ct])
dev.off()
}


#' ### Mean-difference plots
#'
#' Displays log-FCs from the linear model fit against the average log-CPM values.
#+
plotMD(cfitbayes,column=1)
plotMD(cfitbayes,column=2)
plotMD(cfitbayes,column=3)

#' ### Toptables {.tabset .tabset-pills}
#'
#' The following tables summarize results for each contrast.
#'
#+ include=FALSE
tpl <- sapply(mycontr,FUN=function(x){

  t.res <- res2%>%
    dplyr::filter(CONTR==x)%>%
    dplyr::filter(TOPFLAG)%>%
    arrange(pval,desc(abs(logFC)))%>%
    dplyr::select(PARAM,logFC,pval,pval.adj,UPDOWN) #%>%
    #dplyr::mutate_if(is.numeric,round,digits=4)

  return(t.res)
},simplify = FALSE, USE.NAMES = TRUE)

#+
htmltools::tagList(lapply(mycontr,function(x){DT::datatable(tpl[[x]],rownames=FALSE,caption=x)}))

#+
wb <- createWorkbook()
for(i in 1:length(mycontr)){
  addWorksheet(wb, mycontr[i])
  td <-  res2%>%
    dplyr::filter(CONTR==mycontr[i])%>%
    dplyr::filter(TOPFLAG)%>%
    arrange(pval,desc(abs(logFC)))%>%
    dplyr::select(PARAM,logFC,pval,pval.adj,UPDOWN,desc,chr,start,end,size)
  writeData(wb, sheet = i, td, rowNames=FALSE)
}
addWorksheet(wb, "Contrasts")
writeData(wb, sheet = length(mycontr)+1,data.frame(Coefficient=rownames(linfct),linfct), rowNames=FALSE)
saveWorkbook(wb, file = "../output/genelist_de_gsea_immpheno_rnaseq_updated.xlsx", overwrite = TRUE)



#' ### Shared genes {.tabset}
#+ results="asis"
#+ fig.width = 10,fig.height = 6, results="asis"
#+
venn_list <- list( INFLAMED_EXCLUDED_UP = res2 %>% filter(CONTR == "INFLAMED_EXCLUDED",
                                                          TOPFLAG == TRUE,
                                                          UPDOWN == "UP") %>% pull(PARAM),

                   INFLAMED_EXCLUDED_DOWN = res2 %>% filter(CONTR == "INFLAMED_EXCLUDED",
                                                            TOPFLAG == TRUE,
                                                            UPDOWN == "DOWN") %>% pull(PARAM),

                   INFLAMED_DESERT_UP  = res2 %>% filter(CONTR == "INFLAMED_DESERT",
                                                         TOPFLAG == TRUE,
                                                         UPDOWN == "UP") %>% pull(PARAM),

                   INFLAMED_DESERT_DOWN  =  res2 %>% filter(CONTR == "INFLAMED_DESERT",
                                                            TOPFLAG == TRUE,
                                                            UPDOWN == "DOWN") %>% pull(PARAM),

                   EXCLUDED_DESERT_UP = res2 %>% filter(CONTR == "EXCLUDED_DESERT",
                                                        TOPFLAG == TRUE,
                                                        UPDOWN == "UP") %>% pull(PARAM),

                   EXCLUDED_DESERT_DOWN = res2 %>% filter(CONTR == "EXCLUDED_DESERT",
                                                          TOPFLAG == TRUE,
                                                          UPDOWN == "DOWN") %>% pull(PARAM)
)

to_upset <- UpSetR::fromList(venn_list)
upset_plot <- upset(as.data.frame(to_upset), nsets=length(venn_list), order.by = "freq", nintersects=NA, keep.order=TRUE)
upset_plot
wd=getwd()
svg(paste(wd, "Upset_plot.svg", sep="/"), width = 10, height = 8)
upset_plot
dev.off()

#' ### Shared genes table {.tabset .tabset-pills}
#'
#'

x <- get.venn.partitions(venn_list)
x$..values.. <- lapply(x$..values.. , function(x) paste0(x, collapse = ", "))
x <- data.frame(x [ unique(which(rowSums( x[, c(1:which(colnames(x)== "..set..")-1)] )>1)) ,
                    - which(colnames(x)== "..set..")] )   # subsetting table for only shared

if(sum(x$..count.. ==0) > 0) {
  x = x[ ! x$..count.. ==0 ,]
}

#+
htmltools::tagList(DT::datatable(x,rownames=FALSE,caption="Shared genes", options = list(
  autoWidth = TRUE,
  columnDefs = list(list(width = '600px', targets = "..values..")))))



#' ## Signature analysis {.tabset .tabset-pills}
#'
#+

gs.genes <- lapply(gsets,function(x){
  x[["genes"]]
})
idx <- ids2indices(gs.genes,rownames(cfitbayes))

# camera
res.cam <- sapply(mycontr,
                  function(x){
                    res <- camera(v,index=idx,design=design.x,contrast = linfct[,x])
                    res1 <- res%>%
                      rownames_to_column(var="PARAM")
                    res1$CONTR <- x
                    return(res1)
                  },
                  simplify=FALSE
)
res.cam1 <- do.call(rbind,res.cam)

#' ### GSEA Camera
#+
htmltools::tagList(lapply(mycontr,function(x){DT::datatable(
  res.cam1%>%
    dplyr::filter(CONTR==x)%>%
    dplyr::arrange(PValue), #%>%
    #dplyr::mutate_if(is.numeric,round,digits=4),
  rownames=FALSE,caption=x)}))


#+
wb <- createWorkbook()
for(i in 1:length(mycontr)){
  addWorksheet(wb, mycontr[i])
  td <-  res.cam1%>%
    dplyr::filter(CONTR==mycontr[i])%>%
    dplyr::filter(FDR<alpha.sig)%>%
    arrange(PValue)
  writeData(wb, sheet = i, td, rowNames=FALSE)
}
addWorksheet(wb, "Contrasts")
writeData(wb, sheet = length(mycontr)+1,data.frame(Coefficient=rownames(linfct),linfct), rowNames=FALSE)
saveWorkbook(wb, file = "../output/genesetlist_de_gsea_immpheno_rnaseq_updated.xlsx", overwrite = TRUE)


#' ## Single gene boxplots  {.tabset .tabset-pills}
#'
#' The following figures show boxplots of the pre-post differences for the top `r ntop` parameters for each contrast.
#'
#+ fig.width = 20, fig.height = 4, results="asis"

pheno_col <- c(
  INFLAMED = "#d95f02",
  EXCLUDED = "#1b9e77",
  DESERT = "#7570b3"
)

list_gene_plots <-list()

for (ct in mycontr){
  cat('\n')
  cat(paste("### ",ct,"\n"))

    top.param <- res2%>%
      dplyr::filter(CONTR==ct)%>%
      dplyr::arrange(pval)%>%
      dplyr::slice(seq(ntop))%>%
      ungroup%>%
      pull(PARAM)

  emat.sub <- eset.rpkm[fData(eset.rpkm)$symbol%in%top.param, ]
  dd <- convert_eset_to_long(emat.sub)%>%
    dplyr::select(Value,symbol,TISSUE,CD8IMMPH)%>%
    dplyr::rename(logRPKM=Value)

    gp <- dd%>%
      ggplot(.,aes(x=TISSUE,y=logRPKM,col=CD8IMMPH))+
      geom_boxplot()+
      scale_color_manual(values =pheno_col)+
      facet_wrap(~symbol,scales = "free", ncol=ntop)+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))


    list_gene_plots[[ct]]<- gp

    print(gp)

    cat('\n')

}

for(ct in names(list_gene_plots)){
  wd=getwd()
  svg(paste(wd,"/", ct,"_topGenes.svg", sep=""), width = 24, height = 4)
  print(list_gene_plots[[ct]])
  dev.off()
}


#' ## Signature boxplots  {.tabset .tabset-pills}
#'
#' The following figures show boxplots of the pre-post differences for the top `r ntop` parameters for each contrast.
#'
#+ fig.width = 16,fig.height = 8, results="asis"

list_sig_plots <-list()

for (ct in mycontr){
  cat('\n')
  cat(paste("### ",ct,"\n"))

  top.param <- res.cam1%>%
    dplyr::filter(CONTR==ct)%>%
    dplyr::arrange(PValue)%>%
    pull(PARAM)
  top.param2 <- unique(gsub("_UP|_DN","",top.param))[1:ntop]

  signature_scores <- data.frame(signature_scores)
  signature_scores$Signature <- rownames(signature_scores)
  colnames(signature_scores) <- gsub("^X", "",colnames(signature_scores))

  signature_tbl = tidyr::pivot_longer(signature_scores,
                                      cols = where(is.numeric), names_to = "sig_sampleid",
                                      values_to = "sscore") %>% dplyr::left_join(samples %>%
                                                                                   tibble::rownames_to_column("SampleID"), by = c(sig_sampleid = "SampleID"))
  #top.param2 <- c(top.param2, "Endothelial", "Mast","Fibroblast")

  dd <- signature_tbl%>%
    dplyr::filter(Signature %in%top.param2)

  gp <- dd%>%
    ggplot(.,aes(x=TISSUE,y=sscore,col=CD8IMMPH))+
    geom_boxplot()+
    scale_color_manual(values =pheno_col)+
    facet_wrap(~Signature,scales = "free")+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  list_sig_plots[[ct]]<- gp

  print(gp)
  cat('\n')

}


for(ct in names(list_sig_plots)){
  wd=getwd()
  svg(paste(wd,"/", ct,"_topSigs1.svg", sep=""), width = 15, height = 8)
  print(list_sig_plots[[ct]])
  dev.off()
}



#' ## Signature boxplots Excluded  {.tabset .tabset-pills}
#'
#' The following figures show boxplots of the pre-post differences for the top `r ntop` parameters for each contrast.
#'
#+ fig.width = 15,fig.height = 4, results="asis"

list_sig_plots <-list()

for (ct in mycontr){
  cat('\n')
  cat(paste("### ",ct,"\n"))

  top.param <- res.cam1%>%
    dplyr::filter(CONTR==ct)%>%
    dplyr::arrange(PValue)%>%
    pull(PARAM)
  top.param2 <- unique(gsub("_UP|_DN","",top.param))[1:ntop]

  signature_scores <- data.frame(signature_scores)
  signature_scores$Signature <- rownames(signature_scores)
  colnames(signature_scores) <- gsub("^X", "",colnames(signature_scores))

  signature_tbl = tidyr::pivot_longer(signature_scores,
                                      cols = where(is.numeric), names_to = "sig_sampleid",
                                      values_to = "sscore") %>% dplyr::left_join(samples %>%
                                                                                   tibble::rownames_to_column("SampleID"), by = c(sig_sampleid = "SampleID"))
  #top.param2 <- c(top.param2, "Endothelial", "Mast","Fibroblast")
  top.param2 <- c("Endothelial", "Mast","Fibroblast")

  dd <- signature_tbl%>%
    dplyr::filter(Signature %in%top.param2)

  gp <- dd%>%
    ggplot(.,aes(x=TISSUE,y=sscore,col=CD8IMMPH))+
    geom_boxplot()+
    scale_color_manual(values =pheno_col)+
    facet_wrap(~Signature,scales = "free")+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


  list_sig_plots[[ct]]<- gp

  print(gp)
  cat('\n')

}


for(ct in names(list_sig_plots)){
  wd=getwd()
  svg(paste(wd,"/", ct,"_topSigs.svg", sep=""), width = 12, height = 5)
  print(list_sig_plots[[ct]])
  dev.off()
}


#' ### Shared signatures {.tabset}
#' Signatures were considered significant with a PValue<0.05.
#+ results="asis"
#+ fig.width = 6,fig.height = 4, results="asis"


list_shared_sigs <- list()

for (ct in mycontr){

  list_shared_sigs[[ct]] <- res.cam1%>%
                        dplyr::filter(CONTR==ct,
                                      PValue < 0.05)%>%
                        pull(PARAM)

}



to_upset <- UpSetR::fromList(list_shared_sigs)
upset_plot <- upset(as.data.frame(to_upset), nsets=length(venn_list), order.by = "freq", nintersects=NA, keep.order=TRUE)
upset_plot
wd=getwd()
svg(paste(wd, "Upset_signature_plot.svg", sep="/"), width = 6, height = 4)
upset_plot
dev.off()


#' ### Shared sigs table {.tabset .tabset-pills}
#' Signatures were considered significant with a PValue<0.05.
#'

x <- get.venn.partitions(list_shared_sigs)
x$..values.. <- lapply(x$..values.. , function(x) paste0(x, collapse = ", "))
x <- data.frame(x [ unique(which(rowSums( x[, c(1:which(colnames(x)== "..set..")-1)] )>1)) ,
                    - which(colnames(x)== "..set..")] )   # subsetting table for only shared

if(sum(x$..count.. ==0) > 0) {
  x = x[ ! x$..count.. ==0 ,]
}

#+
htmltools::tagList(DT::datatable(x,rownames=FALSE,caption="Shared signatures", options = list(
  autoWidth = TRUE,
  columnDefs = list(list(width = '600px', targets = "..values..")))))



#' ***
#' **SESSION INFO**
#+ collapse = TRUE
sessionInfo()

