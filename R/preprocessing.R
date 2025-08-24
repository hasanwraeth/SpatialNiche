#' Load a 10X Visium HD dataset for niche detection
#'
#' This function loads a Visium HD output directory and generate a list
#' containing two tibbles: barcodes with all bins and images for plotting.
#'
#' @param data.dir Path to the directory containing the binned_outputs directory
#' @param size Specifies the bin sizes to read in - defaults to 8um bin
#'
#' @return A list with two tibbles: \code{barcodes} and \code{images}
#'
#' @importFrom arrow read_parquet
#' @importFrom dplyr rename_with
#' @importFrom rjson fromJSON
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/outs/directory'
#' list.files(data_dir) # Should show binned_outputs folder
#' VisiumHD=Load_Spatial_Data(data.dir = data_dir)
#' Barcodes=VisiumHD$barcodes
#' }
#'
Load_Spatial_Data<-function(
    data.dir,
    size="008um"
    )
{images <- Create_Images(data.dir) %>% bind_rows()
tissue_positions_data.dir <- list.files(paste0(data.dir, "/binned_outputs/square_",size,"/spatial"), pattern = "tissue_positions*", full.names = TRUE)
tissue_positions_file <- basename(tissue_positions_data.dir)
tissue_positions_df <- arrow::read_parquet(tissue_positions_data.dir)%>%
  dplyr::rename_with(~ c("barcode", "tissue", "row", "col", "imagerow", "imagecol"))
data.dir_scales <- file.data.dir(data.dir, paste0("/binned_outputs/square_",size,"/spatial/scalefactors_json.json"))
scales <- rjson::fromJSON(file = data.dir_scales)
data.dir_clusters <- file.data.dir(data.dir, paste0("/binned_outputs/square_",size,"/analysis_csv/clustering/gene_expression_graphclust/clusters.csv"))
images_mod <- images %>% filter(data.dir == data.dir)
if(file.exists(data.dir_clusters))
{clusters <- read.csv(data.dir_clusters)
data.dir_umap <- file.data.dir(data.dir, paste0("/binned_outputs/square_",size,"/analysis_csv/umap/gene_expression_2_components/projection.csv"))
umap <- read.csv(data.dir_umap)
barcodes <- tissue_positions_df %>% mutate(imagerow_scaled = imagerow *
                                        scales$tissue_lowres_scalef, imagecol_scaled = imagecol *
                                        scales$tissue_lowres_scalef, imagerow_scaled_round = round(imagerow *
                                                                                                     scales$tissue_lowres_scalef), imagecol_scaled_round = round(imagecol *
                                                                                                                                                                   scales$tissue_lowres_scalef), tissue = as.factor(tissue)) %>%
  left_join(clusters, by = c(barcode = "Barcode")) %>%
  left_join(umap, by = c(barcode = "Barcode")) %>%
  mutate(height = images_mod$height,
         width = images_mod$width)
}else{
  clusters<-NA
  umap<-NA
  barcodes <- tissue_positions_df %>% mutate(imagerow_scaled = imagerow *
                                          scales$tissue_lowres_scalef, imagecol_scaled = imagecol *
                                          scales$tissue_lowres_scalef, imagerow_scaled_round = round(imagerow *
                                                                                                       scales$tissue_lowres_scalef), imagecol_scaled_round = round(imagecol *
                                                                                                                                                                     scales$tissue_lowres_scalef), tissue = as.factor(tissue)) %>%
    mutate(height = images_mod$height,
           width = images_mod$width)}
cell_id_path <- paste0(data.dir, "/barcode_mappings.parquet")
cell_id_df <- arrow::read_parquet(cell_id_path)
cell_id_df_column=paste0('square_',size)
barcodes$cellid=cell_id_df$cell_id[match(barcodes$barcode,cell_id_df%>%pull(cell_id_df_column))]
return(list(images=images,barcodes=barcodes))}


#' Create a tibble for plotting image
#'
#' This function creates a tibble: images that can be used to create
#' a graphical object (grob) that represents a raster image.
#'
#' @param data.dir Path to the directory containing the spatial directory
#'
#' @return An images tibble
#'
#' @importFrom grid rasterGrob unit
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/outs/directory'
#' list.files(data_dir) # Should show spatial folder
#' images=Create_Images(data.dir = data_dir)
#' }
#'
Create_Images<-function (data.dir)
{image <- Get_Spatial_Files(data.dir, "tissue_lowres_image")
grobs <- grid::rasterGrob(image, width = grid::unit(1.1, "npc"), height = grid::unit(1.1, "npc"), hjust=0.495, vjust=0.495)
images <- tibble(data.dir = factor(data.dir), grob = list(grobs), height = nrow(image), width = ncol(image))
return(images)}


#' Create paths and read different spatial files found in outs directory
#'
#' This function generates paths to the appropriate directory and creates all
#' the spatial files.
#'
#' @param data.dir Path to the directory containing the spatial directory
#' @param type Specifies the type of file required
#'
#' @return A bitmap image
#'
#' @importFrom readbitmap read.bitmap
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/outs/directory'
#' list.files(data_dir) # Should show spatial folder
#' image=Get_Spatial_Files(data.dir = data_dir, type="tissue_lowres_image")
#' }
#'
Get_Spatial_Files<-function (data.dir, type)
{if (type == "tissue_lowres_image") {
  x <- readbitmap::read.bitmap(paste(PATH, "/spatial/tissue_lowres_image.png", sep = ""))}
  if (type == "tissue_hires_image") {
    x <- readbitmap::read.bitmap(paste(PATH, "/spatial/tissue_hires_image.png", sep = ""))}
  if (type == "aligned_fiducials") {
    x <- readbitmap::read.bitmap(paste(PATH, "/spatial/aligned_fiducials.jpg", sep = ""))}
  if (type == "detected_tissue_image") {
    x <- readbitmap::read.bitmap(paste(PATH, "/spatial/detected_tissue_image.jpg", sep = ""))}
  return(x)}



