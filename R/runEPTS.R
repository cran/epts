#' Run the EPTS Shiny Application
#'
#' Launches a Shiny application that provides an interactive user interface 
#' to run the functions provided by the epts package.
#'
#' @examples
#' if (interactive()) {
#'   runEPTS()
#' }
#'
#' @export
runEPTS <- function() {
  appDir <- system.file("shiny-app", package = "epts")
  if (appDir == "") {
    stop("Could not find Shiny app directory. Try re-installing `epts`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}