#' The inverse %in% operator
#'
#' Matching operator (x not in y) oposed to the %in%-operator (x in y)
#'
#' @param x	 vector
#' @param y	 vector of same type as x
#'
#' @return A logical vector
#'
#' @export
#' @concept utility
#'
#' @examples
#' \dontrun{
#' Notx=x`%!in%`y
#' }
#'

`%!in%` <- function(x,y)!(`%in%`(x,y))


#' Compute cell-cell distance based on the spatial coordinates
#'
#' This function creates a large distance matrix d.spatial by calculating the
#' euclidean distance of each cell/spot from every other cell/spot
#'
#' @param coordinates a dataframe in which each row gives the spatial locations/coordinates of each cell/spot
#' @param interaction.range The maximum interaction/diffusion range of ligands. This hard threshold is used to filter out the connections between spatially distant cells
#' @param ratio The conversion factor when converting spatial coordinates from Pixels or other units to Micrometers (i.e.,Microns).
#'
#' For example, setting `ratio = 0.25` indicates that 1 pixel equals 0.25um in the coordinates.
#' For 10X Visium HD, it is the ratio of the theoretical spot size (i.e., 8um) over the number of pixels that span the diameter of a theoretical spot size in the full-resolution image (i.e., 'spot_diameter_fullres' in the 'scalefactors_json.json' file).
#' @param tol The tolerance factor to increase the robustness when comparing the center-to-center distance against the `interaction.range`. This can be the half value of cell/spot size in the unit of um.
#'
#' For example, for 10X visium, `tol` can be set as `8/2`; for slide-seq, `tol` can be set as `10/2`.
#' If the cell/spot size is not known, we provide a function `computeCellDistance` to compute the center-to-center distance. `tol` can be the the half value of the minimum center-to-center distance.
#'
#' @return an object of class "dist" giving the pairwise cell-cell distance
#'
#' @importFrom stats dist
#'
#' @export
#' @concept utility
#'
#' @examples
#' \dontrun{
#' Distance <- Compute_Cell_Distance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
#' }
#'

Compute_Cell_Distance <- function(coordinates, interaction.range = NULL, ratio = NULL, tol = NULL){
  if (ncol(coordinates) == 2) {
    colnames(coordinates) <- c("x_cent","y_cent")
  } else {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }

  d.spatial <- stats::dist(coordinates)
  if (!is.null(ratio)) {
    d.spatial <- d.spatial*ratio
  }

  if(!is.null(interaction.range) & !is.null(tol)){
    message("\n Apply a predefined spatial distance threshold based on the interaction length...")
    d.spatial[d.spatial > (interaction.range + tol)] <- NaN
  }
  return(d.spatial)
}


