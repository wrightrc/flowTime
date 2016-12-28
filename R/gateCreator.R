#' Create the gate environment
#'
#' This environment will allow the gates to be accessible throughout the package and prevent
#' interference with the global environment.
#' @export
gateEnv <- new.env()
load(system.file('extdata', 'Gates', 'defaultGates.Rdata', package = 'flowTime'), envir = gateEnv)

#' A gate for the set of all yeast cells
#'
#' Typically set in FSC.A by SSC.A space to exclued any debris
#'     FSC.A  SSC.A
#'    400000  10000
#'     50000  15000
#'     40000  80000
#'   1250000 600000
#'   2000000 500000
#' @format formal class polygonGate
#' @export
yeastGate <- base::get(x = 'yeastGate', envir = gateEnv)

#' A gate for the set of all diploid singlet yeast cells
#'
#' Typically set in FSC.A by FSC.H space
#' Diploids are typically 5um x 6um ellipsoids while haploids are typically
#' 4um x 4um spheroids. As a result, diploids are longer and you get a
#' larger 'area/volume' FSC.A
#'     FSC.A   FSC.H
#'     80000   80000
#'     70000  100000
#'    800000 1000000
#'   1250000 1250000
#'   1000000 1000000
#' @format formal class polygonGate
#' @export
dipsingletGate <- base::get(x = 'hapsingletGate', envir = gateEnv)

#' A gate for the set of all diploid doublets
#'
#' @format formal class polygonGate
#' @export
dipdoubletGate <- base::get(x = 'hapdoubletGate', envir = gateEnv)

#' A gate for the set of all haploid singlets
#' @format formal class polygonGate
#' @export
hapsingletGate <- base::get(x = 'hapsingletGate', envir = gateEnv)

#' A gate for the set of all haploid doublets
#' @format formal class polygonGate
#' @export
hapdoubletGate <- base::get(x = 'hapdoubletGate', envir = gateEnv)



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
#' @param fileName name of the .Rdata file you would like to save these gates
#' within
#'
#' @return a .RData file in the "extdata" folder of the package containing the
#' specified gates
#' @export
#'
#' @examples
#' loadGates("defaultGates.RData")
#' saveGates()
saveGates <- function(yeastGate = 'yeastGate',
                      dipsingletGate = 'dipsingletGate',
                      dipdoubletGate = 'dipdoubletGate',
                      hapsingletGate = 'hapsingletGate',
                      hapdoubletGate = 'hapdoubletGate',
                      fileName = "defaultGates.RData") {
  if(!exists(c("yeastGate", "hapsingletGate", "hapdoubletGate",
               "dipsingletGate", "dipdoubletGate"))) loadGates()
    save(yeastGate, dipsingletGate, dipdoubletGate, hapsingletGate,
        hapdoubletGate, file = system.file("extdata", "Gates", fileName,
                                        package = "flowTime"))
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
#' @param channels vector containing the channels matching the x and y
#' coordinates above
#'
#' @return a polygon gate object
#' @export
#'
#' @examples
#' polygate(x = c(1,1,10000,10000), y = c(1,10000, 10000, 1), )
polygate <- function(x, y, filterID = "newGate", channels = c("FSC.A",
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

#' Set default gates
#' @description Sets specified file to default gates for use in other
#' functions
#'
#' @param gatesFile the full name of the gates file you would like to set as
#' default (e.g. 'C6Gates.RData')
#'
#' @return overwrites the defaultGates.RData file with the specified gates
#' file
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
#' @description Loads a set of yeast gates into active memory to be used in
#' analysis functions
#' @param gatesFile the gates file to be loaded into memory
#'
#' @return gate objects created in the current environment
#' @export
#'
#' @examples
#' loadGates()
loadGates <- function(gatesFile = "defaultGates.RData") {
  load(system.file("extdata", "Gates", gatesFile, package = "flowTime"),
       envir = gateEnv)
}
