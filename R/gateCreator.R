#' Save a yeast gate set
#'
#' @param yeastGate
#' @param dipsingletGate
#' @param dipdoubletGate
#' @param hapsingletGate
#' @param hapdoubletGate
#' @param file
#'
#' @return
#' @export
#'
#' @examples
saveGates <- function(yeastgate = yeastGate, dipsingletgate = dipsingletGate, dipdoubletgate = dipdoubletGate,
                      hapsingletgate = hapsingletGate, hapdoubletgate = hapdoubletGate, fileName = "defaultGates.Rdata") {
  save(yeastGate, dipsingletGate, dipdoubletGate, hapsingletGate, hapdoubletGate, file = paste0(system.file("extdata/Gates/",
                                                                                                            package = "flowTime"), fileName))
}

#' Create a polygon gate
#'
#' @param x
#' @param y
#' @param filterID
#' @param channels
#'
#' @return
#' @export
#'
#' @examples
polygate <- function(x, y, filterID = "newGate", channels = c("FSC.A", "FSC.H")) {
  if (length(x) != length(y) | !is.numeric(x) | !is.numeric(y)) {
    stop("x coordinate vector must be same length as y coordinate vector")
  }

  gate <- polygonGate(filterId = filterID, .gate = matrix(c(x, y), ncol = 2, nrow = length(x),
                                                          dimnames = list(rep(NA, 5), channels)))
  return(gate)
}

#' Ploidy: Tries to guess the ploidy of a given flowframe using FSC.A/FSC.H ratio.
#' Diploids are typically 5um x 6um ellipsoids while haploids are typically 4um x 4um spheroids. As a result, diploids are longer and you get a larger 'area/volume' FSC.A. 'Width' might also be useful.
#'
#' @param flowframe
#'
#' @return
#' @export
#'
#' @examples
ploidy <- function(flowframe) {
  # Find FSC.A/FSC.H.  This is close to 1 for diploids and close to .8 for haploids Test this
  # assumption!!!!!
  fsca <- summary(flowframe)[4, 1]
  fsch <- summary(flowframe)[4, 7]
  quotient <- fsca/fsch
  if (quotient > 0.92) {
    return(c("Diploid", quotient))
  } else {
    return(c("Haploid", quotient))
  }
}

#'
#' Sets specified file to default gates for use in other functions
#'
#' @param gatesFile
#'
#' @return
#' @export
#'
#' @examples
setGates <- function(gatesFile) {
  load(paste0(system.file("extdata/Gates/", package = "flowTime"), gatesFile))
  saveGates(fileName = paste0(system.file("extdata/Gates/", package = "flowTime"), "defaultGates.RData"))
}

#' Load a yeast gate file
#' into active memory
#' @param gatesFile
#'
#' @return
#' @export
#'
#' @examples
loadGates <- function(gatesFile = "defaultGates.RData") {
  load(paste0(system.file("extdata/Gates/", package = "flowTime"), gatesFile), envir = globalenv())
}
