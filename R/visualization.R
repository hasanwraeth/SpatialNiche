#' Visualize 'features' i.e. genes on a spatial plot
#'
#' This function colors bins on spatial plot according to a 'feature'
#' (i.e. gene expression, PC scores, number of genes detected, etc.)
#'
#' @param barcodes The barcodes tibble
#' @param gene The gene whose expression is to be plotted
#' @param pt.size The dize of the points
#' @param shape The shape of the slice involved (either "circle" or "square")
#'
#' @return A list of ggplot objects
#'
#' @importFrom dplyr pull
#' @importFrom ggplot2 ggplot geom_point coord_cartesian
#' @importFrom scattermore geom_scattermore
#'
#' @export
#' @concept visualization
#'
#' @examples
#' \dontrun{
#' PlotA=Spatial_Feature_Plot(barcodes=Barcodes, gene="GAPDH", pt.size=2,
#'                             shape="circle")
#' }
#'

Spatial_Feature_Plot<-function(barcodes,gene,pt.size=2,shape="circle")
{
  barcodes$Expression<-barcodes %>% pull(gene)

  if(shape=="circle")
  {
    Plot<-barcodes %>%
      ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=Expression)) +
      geom_scattermore(pointsize = pt.size,pixels = rep(2000,2))+
      coord_cartesian(expand = FALSE) +
      xlab("") +
      ylab("") +
      theme_set(theme_bw(base_size = 10))+
      theme_minimal() +
      theme(axis.text = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())+
      scale_color_gradient(low="lightgray",high = "red")+
      labs(color=paste0(gene))

  }else if(shape=="square")
  {
    Plot<-barcodes %>%
      ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,fill=Expression)) +
      geom_point(shape=22,size=pt.size,color=alpha("black",0),stroke=0.25)+
      coord_cartesian(expand = FALSE) +
      xlab("") +
      ylab("") +
      theme_set(theme_bw(base_size = 10))+
      theme_minimal() +
      theme(axis.text = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())+
      scale_fill_gradient(low="lightgray",high = "red")+
      labs(fill=paste0(gene))

  }else{
    stop("Wrong Shape")
  }
  return(Plot)}
