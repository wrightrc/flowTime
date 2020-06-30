#' Save a yeast gate set
#'
#' @param yeastGate a gate object defining the population of yeast cells
#' @param dipsingletGate a gate object defining the population of diploid
#' singlet cells
#' @param dipdoubletGate a gate object defining the population of diploid
#' doublet cells
#' @param hapsingletGate a gate object defining the population of haploid
#' singlet cells
#' @param hapdoubletGate a gate object defining the population of haploid
#' doublet cells
#' @param fileName name of the .Rdata file you would like to save these
#' gates within
#' @param path path to the folder in which you would like to save the gates
#'
#' @return a .RData file in the "extdata" folder of the package containing
#' the specified gates
#' @export
#'
#' @examples
#' loadGates(system.file("extdata/SORPGates.RData", package = "flowTime"))
#' saveGates()
saveGates <- function(yeastGate = NULL,
                      dipsingletGate = NULL,
                      dipdoubletGate = NULL,
                      hapsingletGate = NULL,
                      hapdoubletGate = NULL,
                      path = getwd(),
                      fileName = "defaultGates.RData") {
    save(yeastGate, dipsingletGate, dipdoubletGate, hapsingletGate,
        hapdoubletGate, file = paste(path, fileName, sep = "/"))
}


#' Create a polygon gate
#'
#' @param x a vector of x coordinates
#' @param y a vector of y coordinates
#' @param filterID name of the gate
#' @param channels vector containing the channels matching the x and y
#' coordinates above
#'
#' @return a polygon gate object
#' @export
#'
#' @examples
#' polyGate(x = c(1,1,10000,10000), y = c(1,10000, 10000, 1), )
polyGate <- function(x, y, filterID = "newGate", channels = c("FSC.A",
                                                                "FSC.H")) {
  if (length(x) != length(y) | !is.numeric(x) | !is.numeric(y)) {
    stop("x coordinate vector must be same length as y coordinate vector")
  }

  gate <- polygonGate(filterId = filterID, .gate = matrix(c(x, y), ncol = 2,
                                                    nrow = length(x),
                                                    dimnames = list(rep(NA,
                                                                  length(x)),
                                                                  channels)))
  return(gate)
}

#' Guess the ploidy of a given flowframe
#' @description Use the FSC.A/FSC.H ratio. Diploids are typically 5um x 6um
#' ellipsoids while haploids are typically 4um x 4um spheroids. As a result,
#' diploids are longer and you get a larger 'area/volume' FSC.A. 'Width' might
#' also be useful on certain cytometers.
#'
#' @param flowframe the flowFrame you would like to identify the ploidy of
#'
#' @return "Diploid" or "Haploid" and the mean FSC.A/FSC.H quotient
#' @export
#'
#' @examples dat <- read.flowSet(path = system.file("extdata", "ss_example",
#' package = "flowTime"), alter.names = TRUE)
#' ploidy(dat$A01.fcs)
#'
ploidy <- function(flowframe) {
  # Find FSC.A/FSC.H.  This is close to 1 for diploids and close to .8 for
  # haploids
  # Test this assumption!!!!!
  flowsum <- flowCore::summary(flowframe)
  quotient <- flowsum[4, 1]/flowsum[4, 7]
  if (quotient > 0.92) {
    return(c("Diploid", quotient))
  } else {
    return(c("Haploid", quotient))
  }
}

#' Load a yeast gate file
#' @description Loads a set of yeast gates into active memory to be used in
#' analysis functions
#' @param gatesFile the gates file to be loaded into memory, or path to the
#' gates file
#' @param path The path to the gates file. If 'NULL' this will look through
#' lazy loaded data for the gatesFile
#' @param envir The environment in which to load the gates
#'
#' @return gate objects created in the current environment
#' @export
#' @importFrom utils data
#'
#' @examples
#' loadGates(system.file("extdata/SORPGates.RData", package = "flowTime"))
loadGates <- function(gatesFile = NULL, path = NULL, envir = environment()) {
  if(is.null(path)) data(list = c("dipdoubletGate", "dipsingletGate", "hapdoubletGate",
                                  "hapsingletGate", "yeastGate"), envir = envir)
  else
    load(file = paste(path, gatesFile, sep = "/"), envir = globalenv())
}

