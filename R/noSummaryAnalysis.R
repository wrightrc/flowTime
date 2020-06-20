#' Analysis of steady state fluorescence flow cytometry
#' @description Generates a data frame which can be used to visualize and
#' analyze steady state flow cytometry data. Steady state in this case means
#' that
#' @param flowset your flowSet to be analyzed
#' @param ploidy \code{character} gate to subset your flowset based on the
#' ploidy of you strains
#' @param only \code{character} which population of events to analyze,
#' 'yeast', 'singlets', or 'doublets'?
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
steadyState <- function(flowset, gated = FALSE, ploidy = NA, only = NA) {
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
      message("No ploidy defined, only yeast gate will be applied")
      yeast <- Subset(flowset, yeastGate)
    }
    ### Convert flowSet to dataframe containing all events for each subset
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

#' Generate a tidy dataset from time-course flow cytometry data
#' @description Generates a tibble containing all parameters and phenoData
#' from a flowSet which can be used to visualize and
#' analyze steady state flow cytometry data. Steady state in this case means
#' that
#' @param flowset your flowSet to be analyzed
#' @param ploidy \code{character} gate to subset your flowset based on the
#' ploidy of you strains
#' @param only \code{character} which population of events to analyze,
#' 'yeast', singlets', or 'doublets'?
#' @param gated \code{boolean} is the data already gated?
#' @return a data frame containing all of the selected subset of events from
#' the original flowSet for all parameters including experiment time, etime
#' @export
#' @examples
#' example
#' plate1<-read.flowSet(path=system.file("extdata", "tc_example", package = "flowTime"), alter.names = TRUE)
#' annotation <- read.csv(system.file("extdata", "tc_example.csv", package = "flowTime"))
#' plate1 <- annotateFlowSet(plate1, annotation)
#' tidy_dat <- tidyFlow(plate1, gated = TRUE)
#' head(tidy_dat)
tidyFlow <- function(flowset, gated = FALSE, ploidy = NA, only = NA) {
  tidy_dat <- steadyState(flowset, gated, ploidy, only)
  tidy_dat <- dplyr::rename(tidy_dat, name = X)

  #Generate time columns
  time <- fsApply(flowset, function(frame) {
    btime <- as.numeric(unlist(strsplit(keyword(frame)$`$BTIM`, split = ":")))
    btime <- btime[1] * 60 + btime[2] + btime[3]/60 + btime[4]/6000
    atime <- as.numeric(keyword(frame)$`#ACQUISITIONTIMEMILLI`)/1000/60
    tstep <- as.numeric(keyword(frame)$`$TIMESTEP`)
    name <- gsub(pattern = ".fcs", replacement = "", keyword(frame)$GUID)
    vol <- as.numeric(keyword(frame)$`$VOL`)/1000
    events <- as.numeric(keyword(frame)$`$TOT`)
    return(c(name = name, btime = btime, atime = atime, tstep = tstep,
             vol = vol, events = events))
  })
  #Check function for numeric character columns
  numericcharacters <- function(x) {
    !any(is.na(suppressWarnings(as.numeric(x)))) & is.character(x)
  }
  #Convert numeric characters into numerics
  time <- tibble::as_tibble(time) %>% dplyr::mutate_if(numericcharacters,as.numeric)
  #Join time with tidy_dat
  tidy_dat <- dplyr::left_join(tidy_dat, time)
  #Calculate experiment time
  tidy_dat <- tidy_dat %>% dplyr::mutate(etime = btime - min(btime) + Time*tstep/60)
  tidy_dat
}


