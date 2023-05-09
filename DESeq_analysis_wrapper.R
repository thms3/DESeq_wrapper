#!/usr/bin/env Rscript
#@uthor : Thomas NEFF
#@description : Wrapper for DESeq RNAseq analysis

# lib ---------------------------------------------------------------------
library(DESeq2)
library(xlsx)

# fun ---------------------------------------------------------------------

#' Create DESeq dataset
#' 
#' @param counts reads counts file
#' @param metadata metadata file associated with reads counts file
#' @param var variables for regression
#' @param design formula for design model matrix
#' @param sizefactors factors for normalization
#' 
#' @return  DESeq dataset
dds_dataset <- function(counts,metadata,var,design,sizefactors=NULL){

  if (all(var %in% colnames(metadata))) {
    coldata <- metadata[,var]
    coldata <- sapply(colnames(coldata), function(x) as.factor(coldata[[x]]))
    coldata <- as.data.frame(x = unclass(coldata),stringsAsFactors=T)
  } else {
    stop("Error in names variables")
  }
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,colData = coldata,design=eval(design))
  
  if (is.null(sizefactors)) {
    dds <- DESeq2::estimateSizeFactors(dds)
  } else {
    DESeq2::sizeFactors(dds) <- sizefactors
  }
  
  return(dds)
}

#' Result of DESeq analysis
#' 
#' @param dds DESeq dataset
#' @param contrast var of interest for test DE
#' @param method_padj method for p value adjustement
#' @param no.na delete na resulting after cutoff
#' 
#' @return Large list data=DESeqDataset, res=DESeqResults, res.data.frame= result dataframe
dds_res <- function(dds,contrast,method_padj="BH"){

  res <- DESeq2::results(object = dds,contrast = contrast,pAdjustMethod = method_padj)
  
  res.data.frame <- as.data.frame(res)
  res.data.frame <- res.data.frame[!is.na(res.data.frame$padj),]
  res.data.frame <- res.data.frame[order(res.data.frame[,"padj"]),]
  
  return(list("data"=dds,"res"=res,"res.data.frame"=res.data.frame))
}

#' Save the results in xlsx file
#' 
#' @param dds_lrt result of dds_res() function associated with DESeq2::DESeq(object = dds,test = "LRT",...)
#' @param dds_wald result of dds_res() function associated with DESeq2::DESeq(object = dds, test = "Wald") (default param.)
#' @param file file pathway to save result in xslx format
dds_res_xlsx <- function(dds_wald=NULL,dds_lrt=NULL,file){
  
  #check dds_res() returns 
  stopifnot(!is.null(dds_wald) & is.data.frame(dds_wald$res.data.frame),!is.null(dds_lrt) & is.data.frame(dds_lrt$res.data.frame))
  
  if (file.exists(file)) {
    system(paste0('rm ',file))
  }
  
  if ((!is.null(dds_lrt)) & (!is.null(dds_wald))) {
    write.xlsx(x = dds_lrt$res.data.frame, sheetName = "LRT",file = file)
    write.xlsx(x = dds_wald$res.data.frame, sheetName = "WALD", file = file, append = TRUE)
  } else {
    
    if (!is.null(dds_lrt)) {
      write.xlsx(x = dds_lrt$res.data.frame, sheetName = "LRT",file = file)
    }
    
    if (!is.null(dds_wald)) {
      write.xlsx(x = dds_wald$res.data.frame, sheetName = "WALD",file = file)
    }
    
    if (is.null(dds_lrt) & is.null(dds_wald)) {
      stop("dds_lrt=NULL and dds_wald=NULL")
    }
  }
}
