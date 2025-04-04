#' Launch the Package Shiny Dashboard
#'
#' This function launches the Shiny dashboard included in the package.
#'
#' @export
#' @importFrom shiny runApp
#' @return A Shiny application object
#' @examples
#' if(interactive()) {
#'   launch_dashboard()
#' }
launch_dashboard <- function() {
  app_dir <- system.file("shiny", "results_display", package = "TECAT")
  if (app_dir == "") {
    stop("Could not find app directory. Try reinstalling the package.")
  }
  shiny::runApp(app_dir, display.mode = "normal")
}