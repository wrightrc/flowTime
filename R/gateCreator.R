#' Save a yeast gate set
#'
#' @param yeastgate a gate object defining the population of yeast cells
#' @param dipsingletgate a gate object defining the population of diploid singlet cells
#' @param dipdoubletgate a gate object defining the population of diploid doublet cells
#' @param hapsingletgate a gate object defining the population of haploid singlet cells
#' @param hapdoubletgate a gate object defining the population of haploid doublet cells
#' @param fileName name of the .Rdata file you would like to save these gates within
#'
#' @return a .RData file in the "extdata" folder of the package containing the specified gates
#' @export
#'
#' @examples loadGates()
#' saveGates()
saveGates <- function(yeastgate = yeastGate, dipsingletgate = dipsingletGate, dipdoubletgate = dipdoubletGate,
                      hapsingletgate = hapsingletGate, hapdoubletgate = hapdoubletGate, fileName = "defaultGates.Rdata") {
  save(yeastGate, dipsingletGate, dipdoubletGate, hapsingletGate, hapdoubletGate, file = system.file("extdata", "Gates", fileName, package = "flowTime"))
}

#' List existing gate sets saved in package memory
#'
#' @return lists the gates files saved within the flowTime package
#' @export
#'
#' @examples listGates()
listGates <- function(){
  list.files(system.file("extdata", "Gates", package = "flowTime"))
}

#' Create a polygon gate
#'
#' @param x a vector of x coordinates
#' @param y a vector of y coordinates
#' @param filterID name of the gate
#' @param channels vector containing the channels matching the x and y coordinates above
#'
#' @return a polygon gate object
#' @export
#'
#' @examples
#' polygate(x = c(1,1,10000,10000), y = c(1,10000, 10000, 1), )
polygate <- function(x, y, filterID = "newGate", channels = c("FSC.A", "FSC.H")) {
  if (length(x) != length(y) | !is.numeric(x) | !is.numeric(y)) {
    stop("x coordinate vector must be same length as y coordinate vector")
  }

  gate <- polygonGate(filterId = filterID, .gate = matrix(c(x, y), ncol = 2, nrow = length(x), dimnames = list(rep(NA, length(x)), channels)))
  return(gate)
}

#' Guess the ploidy of a given flowframe
#' @description Use the FSC.A/FSC.H ratio. Diploids are typically 5um x 6um ellipsoids while haploids are typically 4um x 4um spheroids. As a result, diploids are longer and you get a larger 'area/volume' FSC.A. 'Width' might also be useful on certain cytometers.
#'
#' @param flowframe the flowFrame you would like to identify the ploidy of
#'
#' @return "Diploid" or "Haploid" and the mean FSC.A/FSC.H quotient
#' @export
#'
#' @examples dat <- read.flowSet(path = system.file("extdata", "ss_example", package = "flowTime"), alter.names = TRUE)
#' ploidy(dat$A01.fcs)
#'
ploidy <- function(flowframe) {
  # Find FSC.A/FSC.H.  This is close to 1 for diploids and close to .8 for haploids Test this
  # assumption!!!!!
  flowsum <- flowCore::summary(flowframe)
  quotient <- flowsum[4, 1]/flowsum[4, 7]
  if (quotient > 0.92) {
    return(c("Diploid", quotient))
  } else {
    return(c("Haploid", quotient))
  }
}

#' Set default gates
#' @description Sets specified file to default gates for use in other functions
#'
#' @param gatesFile the full name of the gates file you would like to set as default (e.g. 'C6Gates.RData')
#'
#' @return overwrites the defaultGates.RData file with the specified gates file
#' @export
#'
#' @examples
#' loadGates("defaultGates.RData")
#' setGates("defaultGates.RData")
setGates <- function(gatesFile = "defaultGates.RData") {
  load((system.file("extdata", "Gates", gatesFile, package = "flowTime")))
  saveGates(fileName = "defaultGates.RData")
}

#' Load a yeast gate file
#' @description Loads a set of yeast gates into active memory to be used in analysis functions
#' @param gatesFile the gates file to be loaded into memory
#'
#' @return gate objects created in the current environment
#' @export
#'
#' @examples loadGates()
loadGates <- function(gatesFile = "defaultGates.RData") {
  load(system.file("extdata", "Gates", gatesFile, package = "flowTime"), envir = globalenv())
}
