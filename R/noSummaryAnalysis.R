#' Analysis of steady state fluorescence flow cytometry
#' @description Generates a data frame which can be used to visualize and
#' analyze steady state flow cytometry data. Steady state in this case means
#' that
#' @param flowset your flowSet to be analyzed
#' @param ploidy \code{character} gate to subset your flowset based on the
#' ploidy of you strains
#' @param only \code{character} which population of events to analyze,
#' 'yeast', singlets', or 'doublets'?
#' @param gated \code{boolean} is the data already gated?
#' @return a data frame containing all of the selected subset of events from
#' the original flowSet
#' @export
#' @examples
#' dat <- read.flowSet(path = system.file("extdata", "ss_example",
#' package = "flowTime"), alter.names = TRUE)
#' annotation <- read.csv(system.file("extdata", "ss_example.csv",
#' package = "flowTime"))
#' dat <- annotateFlowSet(dat, annotation, mergeBy = "name")
#' loadGates(gatesFile = 'SORPGates')
#' steadyState(dat, gated = FALSE, ploidy = "diploid", only = "singlets")
#'
steadyState <- function(flowset, gated = FALSE, ploidy = "diploid", only = "singlets") {
  ### Number of cells (experiments) in the flowSet
  n_experiment <- length(flowset)

  ### Pulling out data for specific channel to be used
  #channel <- flowSet[,'FL2.A',drop=FALSE]

  ### Gate the samples
  if (gated == FALSE) {
    if (!exists(c("yeastGate", "hapsingletGate", "hapdoubletGate",
                  "dipsingletGate", "dipdoubletGate")))
      loadGates()
    yeastGate <- get("yeastGate")
    hapsingletGate <- get("hapsingletGate")
    hapdoubletGate <- get("hapdoubletGate")
    dipsingletGate <- get("dipsingletGate")
    dipdoubletGate <- get("dipdoubletGate")
    if (ploidy == "haploid") {
      print("Gating with haploid gates...")
      yeast <- Subset(flowset, yeastGate)
      singlets <- Subset(yeast, hapsingletGate)
      doublets <- Subset(yeast, hapdoubletGate)
    } else if (ploidy == "diploid") {
      print("Gating with diploid gates...")
      yeast <- Subset(flowset, yeastGate)
      singlets <- Subset(yeast, dipsingletGate)
      doublets <- Subset(yeast, dipdoubletGate)
    } else {
      stop("Error: You must define ploidy=\"haploid\" or ploidy=\"diploid\"")
    }
    ### Convert flowSet to dataframe containing all events for each subset
    if (only == FALSE) {

      print("Converting all yeast events...")
      yeastdF <- plyr::ddply(pData(yeast), colnames(pData(yeast))[-1],
                             function(tube) {
                               fsApply(x = yeast[tube$name], rbind, use.exprs = TRUE)
                             })
      print("Converting doublets events...")
      doubletsdF <- plyr::ddply(pData(doublets), colnames(pData(doublets))[-1],
                                function(tube) {
                                  fsApply(x = doublets[tube$name], rbind, use.exprs = TRUE)
                                })

      print("Converting singlets events...")
      singletsdF <- plyr::ddply(pData(doublets), colnames(pData(doublets))[-1],
                                function(tube) {
                                  fsApply(x = doublets[tube$name], rbind, use.exprs = TRUE)
                                })
    }
    if (only == "singlets") {
      print("Converting singlets events...")
      singletsdF <- plyr::ddply(pData(singlets), colnames(pData(singlets))[-1],
                                function(tube) {
                                  fsApply(x = singlets[tube$name], rbind, use.exprs = TRUE)
                                })
      return(singletsdF)
    } else if (only == "doublets") {
      print("Converting doublets events...")
      doubletsdF <- plyr::ddply(pData(doublets), colnames(pData(doublets))[-1],
                                function(tube) {
                                  fsApply(x = doublets[tube$name], rbind, use.exprs = TRUE)
                                })
      return(doubletsdF)
    } else if (only == "yeast") {
      print("Converting all yeast events...")
      yeastdF <- plyr::ddply(pData(yeast), colnames(pData(yeast))[-1],
                             function(tube) {
                               fsApply(x = yeast[tube$name], rbind, use.exprs = TRUE)
                             })
      return(yeastdF)
    }
  }
  else {
    print("Converting all events...")
    yeastdF <- plyr::ddply(pData(flowset), colnames(pData(flowset))[-1],
                           function(tube) {
                             fsApply(x = flowset[tube$name], rbind, use.exprs = TRUE)
                           })
    return(yeastdF)
  }
}


