#' Create a character vector of a square region for sub-setting
#'
#' This function uses a given size (in microns) whose center is a given barcode
#' to create a character vector of barcodes from the barcodes tibble that can be
#' used to subset a square region of the Visium HD dataset.
#'
#' @param spot Barcode at the center of the square to be generated
#' @param size_microns The length of the square in microns
#' @param barcodes The barcodes tibble
#' @param binsize The binsize of the current analysis (default=8)
#'
#' @return A character vector of the barcodes in the square region
#'
#' @importFrom dplyr pull filter
#'
#' @export
#' @concept slicing
#'
#' @examples
#' \dontrun{
#' SquareA=Get_Square(spot='s_008um_00128_00230-1', size_microns=350,
#'                     barcodes=Barcodes, binsize=8)
#' }
#'

Get_Square<-function(spot,size_microns,barcodes,binsize=8)
{Xcenter<-barcodes$col[match(spot,barcodes$barcode)]
  Ycenter<-barcodes$row[match(spot,barcodes$barcode)]
  AddFactor<-round(size_microns/(2*binsize))
  Xmin<-Xcenter-AddFactor
  Xmax<-Xcenter+AddFactor
  Ymin<-Ycenter-AddFactor
  Ymax<-Ycenter+AddFactor
  SquareSection<-barcodes %>% dplyr::filter(col >= Xmin & col <= Xmax & row >= Ymin & row <= Ymax) %>% dplyr::pull(barcode)
  return(SquareSection)}

#' Create a character vector of a rectangular region for sub-setting
#'
#' This function uses a given size (in microns) whose center is a given barcode
#' to create a character vector of barcodes from the barcodes tibble that can be
#' used to subset a rectangular region of the Visium HD dataset.
#'
#' @param spot Barcode at the center of the rectangle to be generated
#' @param size_microns The length of the rectangle in microns
#' @param xfactor The factor of the horizontal line of the rectangle
#' @param yfactor The factor of the vertical line of the rectangle
#' @param barcodes The barcodes tibble
#' @param binsize The binsize of the current analysis (default=8)
#'
#' @return A character vector of the barcodes in the rectangular region
#'
#' @importFrom dplyr pull filter
#'
#' @export
#' @concept slicing
#'
#' @examples
#' \dontrun{
#' RectangleA=Get_Rectangle(spot='s_008um_00128_00230-1', size_microns=350,
#'                     barcodes=Barcodes, binsize=8)
#' }
#'

Get_Rectangle<-function(spot,size_microns,xfactor,yfactor,barcodes,binsize=8)
{Xcenter<-barcodes$col[match(spot,barcodes$barcode)]
  Ycenter<-barcodes$row[match(spot,barcodes$barcode)]
  AddFactor<-round(size_microns/(2*binsize))
  Xmin<-Xcenter-(AddFactor*xfactor)
  Xmax<-Xcenter+(AddFactor*xfactor)
  Ymin<-Ycenter-(AddFactor*yfactor)
  Ymax<-Ycenter+(AddFactor*yfactor)
  RectSection<-barcodes %>% dplyr::filter(col >= Xmin & col <= Xmax & row >= Ymin & row <= Ymax) %>%  dplyr::pull(barcode)
  return(RectSection)}

#' Create a character vector of a circular region for sub-setting
#'
#' This function uses a given radius size (in microns) whose center is a given
#' barcode to create a character vector of barcodes from the barcodes tibble
#' that can be used to subset a circular region of the Visium HD dataset.
#'
#' @param spot Barcode at the center of the circle to be generated
#' @param size_microns The radius of the circle in microns
#' @param barcodes The barcodes tibble
#' @param data.dir Path to the directory containing the binned_outputs directory
#' @param celltype The type of cell at the center
#' @param binsize The binsize of the current analysis (default="008um")
#'
#' @return A character vector of the barcodes in the circular region
#'
#' @export
#' @concept slicing
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/outs/directory'
#' list.files(data_dir) # Should show binned_outputs folder
#' CircleA=Get_Circle(spot='s_008um_00128_00230-1', size_microns=350,
#'                     barcodes=Barcodes, data.dir=data_dir, celltype=NA,
#'                     binsize="008um")
#' }
#'

Get_Circle<-function(spot,size_microns,barcodes,data.dir,celltype=NA,binsize="008um")
{data.dir_scales <- paste0(data.dir, "/binned_outputs/square_",binsize,"/spatial/scalefactors_json.json")
  scales <- rjson::fromJSON(file = data.dir_scales)
  Scale<-(size_microns*scales$spot_diameter_fullres)/as.numeric(unlist(strsplit(binsize,"um"))[1])
  Index<-match(spot,barcodes$barcode)
  Result<-vector("list",length=length(Index))
  for(jj in 1:length(Index))
  {Distance<-sqrt(((barcodes$imagecol-barcodes$imagecol[Index[jj]])^2) + ((barcodes$imagerow-barcodes$imagerow[Index[jj]])^2))
  barcodes$Distance<-Distance
  if(!is.na(celltype))
  {ValTh <- sum(barcodes$DeconvolutionLabel1[barcodes$Distance<min(Scale)]==celltype,na.rm = T)
  if(ValTh < 25)
  {next}}
  if(length(Scale)>1)
  {Result[[jj]]<-lapply(Scale,function(X){return(barcodes$barcode[barcodes$Distance < X])})
  }else{
    Result[[jj]]<-barcodes$barcode[barcodes$Distance < Scale]}}
  if(length(Scale)>1)
  {Rxx<-vector("list",length=length(Scale))
  names(Rxx)<-as.character(size_microns)
  for(ii in 1:length(Scale))
  {Rxx[[ii]]<-lapply(Result,function(X){return(X[[ii]])})
  Rxx[[ii]]<-unique(unlist(Rxx[[ii]]))}
  return(Rxx)
  }else{
    Result<-unique(unlist(Result))
    return(Result)}}

#' Create a character vector of peripheral cells of a cell based niche for
#' sub-setting
#'
#' This function uses a given size (in microns) whose center is a given cell
#' type to create a character vector of barcodes from the barcodes tibble that
#' can be used to subset peripheral cells of a niche from the Visium HD dataset.
#'
#' @param spot Barcode at the center of the rectangle to be generated
#' @param size_microns The radius of the circle in microns
#' @param barcodes The barcodes tibble
#' @param data.dir Path to the directory containing the binned_outputs directory
#' @param celltype The type of cell at the center
#' @param binsize The binsize of the current analysis (default="008um")
#'
#' @return A character vector of the barcodes surrounding a cell type
#'
#' @importFrom dplyr filter
#'
#' @export
#' @concept slicing
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/outs/directory'
#' list.files(data_dir) # Should show binned_outputs folder
#' PeripheryA=Get_Periphery(barcodes=Barcodes, celltype=Cell_Type, distance= 50,
#'                     data.dir = data_dir, binsize="008um")
#' }
#'

Get_Periphery<-function(barcodes,celltype,distance=50,data.dir, binsize="008um")
{SelectedBarcodes<-barcodes %>% dplyr::filter(DeconvolutionLabel1==celltype)
  Result<-Get_Circle(SelectedBarcodes$barcode,distance,barcodes,data.dir,celltype=celltype, binsize="008um")
  if(length(distance)>1)
  {for(jj in 1:length(Result))
  {Result[[jj]]<-Result[[jj]][Result[[jj]]%!in%SelectedBarcodes$barcode]}
    return(Result)
  }else{
    Result<-Result[Result%!in%SelectedBarcodes$barcode]
    return(Result)}}

#' Create a character vector of peripheral cells of a gene based niche for
#' sub-setting
#'
#' This function uses a given size (in microns) whose center is a cell
#' expressing a particualr gene to create a character vector of barcodes from
#' the barcodes tibble that can be used to subset peripheral cells of a niche
#' from the Visium HD dataset.
#'
#' @param spot Barcode at the center of the rectangle to be generated
#' @param size_microns The radius of the circle in microns
#' @param barcodes The barcodes tibble
#' @param data.dir Path to the directory containing the binned_outputs directory
#' @param gene The gene at the center
#' @param binsize The binsize of the current analysis (default="008um")
#' @param seurat The Seurat object of Visium HD analysis
#'
#' @return A character vector of the barcodes surrounding a cell expressing a
#' gene
#'
#' @importFrom dplyr filter
#'
#' @export
#' @concept slicing
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/outs/directory'
#' list.files(data_dir) # Should show binned_outputs folder
#' PeripheryA=Get_Gene_Periphery(barcodes=Barcodes, gene="GAPDH", distance= 50,
#'                     data.dir = data_dir, binsize="008um", seurat=Seurat)
#' }
#'

Get_Gene_Periphery<-function(barcodes,gene,distance=50,data.dir, binsize="008um",seurat)
{barcodes<-Add_Expression(barcodes,seurat,gene)
  barcodes$DeconvolutionLabelg=barcodes$DeconvolutionLabel1
  barcodes$DeconvolutionLabelg=barcodes%>%pull(gene)>0
  SelectedBarcodes<-barcodes %>% dplyr::filter(DeconvolutionLabelg==TRUE)
  Result<-Get_Circle(SelectedBarcodes$barcode,distance,barcodes,data.dir,binsize="008um")
  if(length(distance)>1)
  {for(jj in 1:length(Result))
  {Result[[jj]]<-Result[[jj]][Result[[jj]]%!in%SelectedBarcodes$barcode]}
    return(Result)
  }else{
    Result<-Result[Result%!in%SelectedBarcodes$barcode]
    return(Result)}}
