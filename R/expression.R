#' Add gene expression data to the barcodes tibble
#'
#' This function uses a Seurat object to load specific gene expression data to
#' the barcodes tibble.
#'
#' @param barcodes The barcodes tibble
#' @param seurat The Seurat object of Visium HD analysis
#' @param genes A list of genes whose expression is to be added
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
#' Barcodes=Add_Expression(barcodes=Barcodes,seurat=Seurat)
#' }
#'

Add_Expression<-function(barcodes,seurat,genes)
{Exp<-SeuratObject::FetchData(seurat,genes)
for(Gx in genes)
{barcodes[,Gx]<-NA
barcodes[match(rownames(Exp),barcodes$barcode),Gx]<-Exp[,Gx]}
return(barcodes)}
