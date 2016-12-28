#' Create an annotation dataframe
#' @description Creates a data frame with rows containing the sample names of
#' your flow set that can then be filled in with experimental metadata.
#'
#' @param yourFlowSet the flowSet to create an annotation data frame for
#'
#' @return annotation_df a data frame containing the sample names of your
#' flow set
#' @export
#'
#' @examples dat <- read.flowSet(path = system.file("extdata", "ss_example",
#' package = "flowTime"), alter.names = TRUE)
#' annotation <- createAnnotation(yourFlowSet = dat)
#' head(annotation)
createAnnotation <- function(yourFlowSet) {
  annotation_df <- data.frame(name = sampleNames(yourFlowSet))
}

#' Annotate a flowSet with experimental metadata
#'
#' @description Add annotations to a flowSets phenoData and plate numbers,
#' strain names, and treatment also set T0
#'
#' @param yourFlowSet a flowSet with sampleNames of the format 'plate#_Well',
#' we typically use the following code chunk to read data from individual
#' plates as exported from BD Accuri C6 software.
#' @param annotation_df A data frame with columns 'well', 'strain',
#' 'treatment', containing all of the wells in the flowset labeled with the
#' strain and treatment in that well.
#' @param mergeBy the unique identifier column
#' @return An annotated flowSet
#' @export
#'
#' @examples dat <- read.flowSet(path = system.file("extdata", "ss_example",
#' package = "flowTime"), alter.names = TRUE)
#' annotation <- read.csv(system.file("extdata", "ss_example.csv", package =
#' "flowTime"))
#' annotateFlowSet(dat, annotation, mergeBy = "name")
#'
annotateFlowSet <- function(yourFlowSet, annotation_df, mergeBy = "name") {
  pData(yourFlowSet)$name <- sampleNames(yourFlowSet)
  mdata <- plyr::join(pData(yourFlowSet), annotation_df, by = mergeBy)
  rownames(mdata) <- as.character(mdata[, mergeBy])
  pData(yourFlowSet) <- mdata
  yourFlowSet
}






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
#' loadGates(gatesFile = 'SORPGates.RData')
#' steadyState(dat, gated = FALSE, ploidy = "diploid", only = "singlets")
#'
steadyState <- function(flowset, gated = FALSE, ploidy = "diploid", only = "singlets") {

  if(!exists(c("yeastGate", "hapsingletGate", "hapdoubletGate",
                "dipsingletGate", "dipdoubletGate"), envir = gateEnv)) loadGates()
  ### Number of cells (experiments) in the flowSet
  n_experiment <- length(flowset)

  ### Pulling out data for specific channel to be used
  #channel <- flowSet[,'FL2.A',drop=FALSE]

  ### Gate the samples
  if (gated == FALSE) {
    if (!exists(c("yeastGate", "hapsingletGate", "hapdoubletGate",
                  "dipsingletGate", "dipdoubletGate"), envir = gateEnv))
      loadGates()
    if (ploidy == "haploid") {
      print("Gating with haploid gates...")
      yeast <- Subset(flowset, flowTime::yeastGate)
      singlets <- Subset(yeast, flowTime::hapsingletGate)
      doublets <- Subset(yeast, flowTime::hapdoubletGate)
    } else if (ploidy == "diploid") {
      print("Gating with diploid gates...")
      yeast <- Subset(flowset, flowTime::yeastGate)
      singlets <- Subset(yeast, flowTime::dipsingletGate)
      doublets <- Subset(yeast, flowTime::dipdoubletGate)
    } else {
      stop("Error: You must define ploidy=\"haploid\" or ploidy=\"diploid\"")
    }
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


