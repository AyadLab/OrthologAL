#' Launch OrthologAL
#' @examples
#' RunOrthologAL()
#' @import shiny
#' @import Seurat
#' @import biomaRt
#' @import ggplot2
#' @import data.table
#' @import viridis
#' @import dplyr
#' @import bslib
#' @import htmlwidgets
#' @import bslib
#' @import shinydashboard
#' @import shinythemes
#' @import DT
#' @export
RunOrthologAL <- function() {
  options(shiny.maxRequestSize = 80000 * 1024^2) #to increase the upload size of seurat object
  shiny::runApp(appDir = system.file('shiny', package = "OrthologAL"))
}
