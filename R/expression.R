#' Add gene expression data to the barcodes tibble
#'
#' This function adds gene expression data from a Seurat object to the
#' previously generated barcodes tibble.
#'
#' @param barcodes The barcodes tibble
#' @param gene The gene at the center
#' @param seurat The Seurat object of Visium HD analysis
#'
#' @return A \code{barcodes} tibble with the gene expression data added to it
#'
#' @importFrom SeuratObject FetchData
#'
#' @export
#' @concept expression
#'
#' @examples
#' \dontrun{
#' Barcodes=Add_Expression(barcodes=Barcodes, gene="GAPDH", seurat=Seurat)
#' }
#'

Add_Expression<-function(barcodes,seurat,gene)
{seurat<-subset(seurat,cells=barcodes$barcode)
  Exp<-FetchData(seurat,gene)
  for(Gx in gene)
  {barcodes[,Gx]<-NA
  barcodes[match(rownames(Exp),barcodes$barcode),Gx]<-Exp[,Gx]}
  return(barcodes)}
