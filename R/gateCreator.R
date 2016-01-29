# #########################
# ###  Cytometer Gates  ###
# #########################
#
#
#
# ###These stricter gates were set for the NemLab cytometer "Special Snowflake" 7sept2015 by EPJ
#
# ### EPJyeastGate
# ### Defines an SSC.A vs FSC.A gate.  Includes only the yeast population from a flowSet
#
# yeastGate <<- polygonGate(filterId="Yeast",
#                              .gate=matrix(c(
#                                #xvalues
#                                4e5,50000,40000,1250000,2000000,
#                                #yvalues
#                                10000,15000,0.8e5,6e5,5e5),
#                                ncol=2,nrow=5,dimnames=list(c("1","1","1","1","1"),c("FSC.A","SSC.A"))))
#
# ### Diploid Gates
# ### Diploids are slightly larger and have better separation between singlets/doublets
#
# dipsingletGate <<- polygonGate(filterId="DipSingletGate",
#                                   .gate=matrix(c(
#                                     #x values
#                                     8e4,7e4,8e5,1.25e6,1e6,
#                                     #y values
#                                     8e4,1e5,1e6,1.25e6,1e6),
#                                     ncol=2,nrow=5,dimnames=list(rep(NA,5),c("FSC.A","FSC.H"))))
#
# dipdoubletGate <<- polygonGate(filterId="DipDoubletGate",
#                                   .gate=matrix(c(
#                                     #x values
#                                     4e5,4.5e4,2e5,1.3e6,2e6,2e6,
#                                     #y values
#                                     2.5e5,5e4,2e5,1.3e6,1.5e6,1e6),
#                                     ncol=2,nrow=6,dimnames=list(rep(NA,6),c("FSC.A","FSC.H"))))
#
# ### Haploid Gates -->>still  need to set (next time I process haploid data)
#
# hapsingletGate <<- polygonGate(filterId="HaploidSingletGate",
#                                   .gate=matrix(c(
#                                     #x values
#                                     1e4,6.5e4,15e5,12e5,4e4,
#                                     #y values
#                                     4e4,8e4,16e5,18e5,4e5),
#                                     ncol=2,nrow=5,dimnames=list(rep(NA,5),c("FSC.A","FSC.H"))))
#
# hapdoubletGate <<- polygonGate(filterId="HaploidDoubletGate",
#                                   .gate=matrix(c(
#                                     #x values
#                                     1e4,2e5,30e5,25e5,25e5,6.75e5,
#                                     #y values
#                                     2e4,5e4,20e5,30e5,5e5,8.5e5),
#                                     ncol=2,nrow=6,dimnames=list(rep(NA,6),c("FSC.A","FSC.H"))))

# save(yeastGate, dipsingletGate, dipdoubletGate, hapsingletGate, hapdoubletGate, file = 'data/EPJGates.RData')
#' Title
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
saveGates <- function(yeastGate = yeastGate, dipsingletGate = dipsingletGate, dipdoubletGate = dipdoubletGate, hapsingletGate = hapsingletGate, hapdoubletGate = hapdoubletGate, file = 'defaultGates.Rdata'){
  save(yeastGate, dipsingletGate, dipdoubletGate, hapsingletGate, hapdoubletGate, file = system.file("extdata", file, package = "flowTime"))
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
polygate <- function(x,y,filterID="newGate",channels=c("FSC.A","FSC.H")) {
  if( length(x) != length(y) | !is.numeric(x) | !is.numeric(y)) {
    stop("x coordinate vector must be same length as y coordinate vector")
  }

  gate <- polygonGate(filterId=filterID,
                      .gate=matrix(c(x,y),
                                   ncol=2,nrow=length(x),dimnames=list(rep(NA,5),channels)))
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
  # Find FSC.A/FSC.H.  This is close to 1 for diploids and close to .8 for haploids
  # Test this assumption!!!!!
  fsca <- summary(flowframe)[4,1]
  fsch <- summary(flowframe)[4,7]
  quotient <- fsca/fsch
  if(quotient>0.92) {
    return(c("Diploid",quotient))
  } else{
    return(c("Haploid",quotient))
  }
}

setGates <- function()
