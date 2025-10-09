#' Load a RCTD result for deconvolution
#'
#' This function loads a sapcexr RCTD result rds object with the appropriate
#' weights that can be used for deconvolution 10X Visium HD dataset.
#'
#' @param data A sapcexr RCTD result rds object
#'
#' @return A list with \code{RCTD_df} tibble and \code{Weights} matrix
#'
#' @export
#' @concept deconvolution
#'
#' @examples
#' \dontrun{
#' data <- 'path/to/RCTD_rds/directory'
#' Deconvolution=Load_RCTD(data = data)
#' }
#'

Load_RCTD<-function(data)
{RCTD<-readRDS(data)
Info=RCTD@results$results_df
weights <- Get_Doublet_Weights_Modified(RCTD@results$results_df, RCTD@results$weights_doublet, RCTD@cell_type_info$info[[2]])
Results <- list(RCTD_df = Info, Weights = t(weights))
return(Results)}

#' Create a Weights matrix from RCTD results
#'
#' This function creates a modified weights matrix using the barcodes,
#' cell_types and doublet_weights.
#'
#' @param results_df The RCTD results tibble
#' @param weights The doublet weights stored in RCTD
#' @param celltypes The cell types of the reference
#'
#' @return A modified \code{Weights} matrix
#'
#' @export
#' @concept deconvolution
#'
#' @examples
#' \dontrun{
#' Weights=Get_Doublet_Weights_Modified(results_df,weights,celltypes)
#' head(Weights)
#' }
#'

Get_Doublet_Weights_Modified <- function(results_df,weights,celltypes) {
  barcodes <- rownames(results_df)
  mod_weights <- matrix(0, nrow = length(barcodes), ncol = length(celltypes))
  rownames(mod_weights) <- barcodes
  colnames(mod_weights) <- celltypes
  indexRow_Certain<-which(results_df$spot_class %in% c('singlet', 'doublet_certain'))
  indexCol_Certain<-match(results_df[indexRow_Certain,'first_type'],colnames(mod_weights))
  mod_weights[cbind(indexRow_Certain,indexCol_Certain)] <- weights[indexRow_Certain,1]
  indexRow_Doublet<-which(results_df$spot_class == "doublet_certain")
  indexCol_Doublet<-match(results_df[indexRow_Doublet,'second_type'],colnames(mod_weights))
  mod_weights[cbind(indexCol_Doublet)] <- weights[indexRow_Doublet,2]
  return(mod_weights)}

#' Add RCTD results and weights to the barcodes tibble
#'
#' This function adds the RCTD deconvolution results to the previously generated
#' barcodes tibble.
#'
#' @param barcodes The barcodes tibble
#' @param deconv The loaded RCTD deconvolution object
#' @param addWeights Add modified doublet weights to the tibble
#'
#' @return A \code{barcodes} tibble with the deconvolution results added to it
#'
#' @export
#' @concept deconvolution
#'
#' @examples
#' \dontrun{
#' Barcodes=Add_Deconvolution_Info(barcodes=Barcodes,deconv=Deconvolution)
#' }
#'

Add_Deconvolution_Info<-function(barcodes,deconv,addWeights=FALSE)
{ResultsDF<-deconv$RCTD_df
index<- match(rownames(ResultsDF),barcodes$barcode)
barcodes$DeconvolutionClass<-NA
barcodes$DeconvolutionClass[index]<-as.vector(ResultsDF$spot_class)
barcodes$DeconvolutionLabel1<-NA
barcodes$DeconvolutionLabel1[index]<-as.vector(ResultsDF$first_type)
barcodes$DeconvolutionLabel2<-NA
barcodes$DeconvolutionLabel2[index]<-as.vector(ResultsDF$second_type)
if(addWeights)
{Weights<-deconv$Weights
index<- match(colnames(Weights),barcodes$barcode)
Names<-gsub(" ","",rownames(Weights))
for(jj in 1:nrow(Weights))
{barcodes[,Names[jj]]<-NA
barcodes[index,Names[jj]]<-Weights[jj,]}}
return(barcodes)}
