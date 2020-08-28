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
  if (gated == FALSE) {
    if (ploidy == "haploid") {
      if (only == FALSE | only == "yeast") {
        if (exists("yeastGate")){
          print("Gating with haploid yeast gate...")
          subset <- Subset(flowset, yeastGate)
        }
        else {
          print("`yeastGate` object not found in environment. Load a
                   gateSet with loadGates, create a `yeastGate` filter object,
                   or set `gated = FALSE` in your call to `summarizeFlow`")
          stop()
        }
      }
      else if(only == "singlets"){
        if (exists("yeastGate") & exists("hapsingletGate")){
          print("Gating with haploid singlet gates...")
          subset <- Subset(flowset, yeastGate & hapsingletGate)
        }
        else {
          print("`yeastGate`  or `hapsingletGate` object not found in
                  environment. Load a
                  gateSet with loadGates, create a `yeastGate` filter object,
                  or set `gated = FALSE` in your call to `summarizeFlow`")
          stop()
        }
      }
      else if(only == "doublets"){
        if (exists("yeastGate") & exists("hapdoubletGate")){
          print("Gating with haploid doublet gates...")
          subset <- Subset(flowset, hapdoubletGate)
        }
        else {
          print("`yeastGate` or `hapdoubletGate` object not found in
                  environment. Load a
                  gateSet with loadGates, create a `yeastGate` filter object,
                  or set `gated = FALSE` in your call to `summarizeFlow`")
          stop()
        }
      }
      else {
        print("`only` value not identified. No further gating applied.")
        subset <- flowset
      }
    }
    else if (ploidy == "diploid") {
      if (only == FALSE | only == "yeast") {
        if (exists("yeastGate")){
          print("Gating with diploid yeast gate...")
          subset <- Subset(flowset, yeastGate)
        }
        else {
          print("`yeastGate` object not found in environment. Load a
                   gateSet with loadGates, create a `yeastGate` filter object,
                   or set `gated = FALSE` in your call to `summarizeFlow`")
          stop()
        }
      }
      else if(only == "singlets"){
        if (exists("yeastGate") & exists("dipsingletGate")){
          print("Gating with diploid singlet gates...")
          subset <- Subset(flowset, yeastGate & dipsingletGate)
        }
        else {
          print("`yeastGate`  or `dipsingletGate` object not found in
                  environment. Load a
                  gateSet with loadGates, create a `yeastGate` filter object,
                  or set `gated = FALSE` in your call to `summarizeFlow`")
          stop()
        }
      }
      else if(only == "doublets"){
        if (exists("yeastGate") & exists("dipdoubletGate")){
          print("Gating with diploid doublet gates...")
          subset <- Subset(flowset, dipdoubletGate)
        }
        else {
          print("`yeastGate` or `dipdoubletGate` object not found in
                  environment. Load a
                  gateSet with loadGates, create a `yeastGate` filter object,
                  or set `gated = FALSE` in your call to `summarizeFlow`")
          stop()
        }
      }
      else {
        print("`only` value not identified. No further gating applied.")
        subset <- flowset
      }
    }
    else {
      print("No ploidy specified. No further gating applied.")
      subset <- flowset
    }
  }
  else if(gated == TRUE) {
    print("No further gating applied.")
    subset <- flowset
  }
  else{
    warning("Unidentified `gated` value. No further gating applied.")
    subset <- flowset
  }
  print("Converting events...")
    dF <- plyr::ddply(pData(subset), colnames(pData(subset))[-1],
                             function(tube) {
                               fsApply(x = subset[tube$name], rbind, use.exprs = TRUE)
                             })
  return(dF)
}

#' Generate a tidy dataset from time-course flow cytometry data
#' @description Generates a tibble containing all parameters and phenoData
#' from a flowSet which can be used to visualize and
#' analyze timecourse flow cytometry data.
#' @param flowset your flowSet to be analyzed
#' @param ploidy \code{character} gate to subset your flowset based on the
#' ploidy of you strains
#' @param only \code{character} which population of events to analyze,
#' 'yeast', singlets', or 'doublets'?
#' @param gated \code{boolean} is the data already gated?
#' @return a data frame containing all of the selected subset of events from
#' the original flowSet for all parameters including experiment time, etime,
#' the time after the initial reading at which each event was collected.
#' @export
#' @examples
#' plate1<-read.flowSet(path=system.file("extdata", "tc_example",
#' package = "flowTime"), alter.names = TRUE)
#' annotation <- read.csv(system.file("extdata", "tc_example.csv",
#' package = "flowTime"))
#' plate1 <- annotateFlowSet(plate1, annotation)
#' tidy_dat <- tidyFlow(plate1, gated = TRUE)
#' head(tidy_dat)
#'
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


