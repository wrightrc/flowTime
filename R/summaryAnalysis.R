
#' Get summary statistics for fluorescence or other data channels of a flowSet
#'
#' @import flowCore
#'
#' @param flowset the \code{flowSet} to create summary statistics for
#' @param channel option \code{character vector} of the data channel(s) to
#' summarize. By default all channels will be summarized. Setting channels does
#' not reduce computation time
#'
#' @return A \code{data frame} containing summary statistics (mean, median,
#' SD) for the specified fluorescent channel and time moments of the flowSet.
#' @export
#'
#' @examples
#' plate1 <- read.flowSet(path = system.file("extdata",
#' "ss_example", package = "flowTime"), alter.names = TRUE)
#' flsummary(flowset = plate1, channel = "FL1.A")
flsummary <- function(flowset, channel) {
  # Number of cells (experiments) in the flowSet
  n_experiment <- length(flowset)

  # Initialize empty matrices/data frames to increase efficiency
  warnings <- c()

  # Get time of each frame in minutes of the day
  btime_raw <- fsApply(flowset, function(x)
    as.numeric(unlist(strsplit(keyword(x)$`$BTIM`, split = ":"))))
  btime <- if(ncol(btime_raw) == 4) apply(btime_raw, 1, function(x) x[1] * 60 +
                      x[2] + x[3]/60 + x[4]/6000) else
           if(ncol(btime_raw) == 3) apply(btime_raw, 1, function(x) x[1] * 60 +
                      x[2] + x[3]/60) else
           stop("Invalid BTIM paramater in FCS file")

  stopifnot(names(btime) == sampleNames(flowset))
  btime <- unname(btime)
  time <- btime - min(btime)

  # Acquisition time - how long it took to take the sample, in seconds
  if (!is.null(keyword(flowset[[1]])$`#ACQUISITIONTIMEMILLI`)) {
    atime <- fsApply(flowset, function(x) {
      as.numeric(keyword(x)$`#ACQUISITIONTIMEMILLI`)/1000
      })
    } else {
    atime <- fsApply(flowset, function(x) {
      max(exprs(x)[,"Time"])/60000
      })
    }

  events <- fsApply(flowset, function(x) length(x[, 1]), use.exprs = TRUE)
  if (!is.null(keyword(flowset[[1]])$`$VOL`)) {
    uL <- fsApply(flowset, function(x) as.integer(keyword(x)$`$VOL`)/1000)
    conc <- events/uL
  } else
    conc <- rep(NA, length(flowset))


  for (i in 1:n_experiment) {
    if (events[i] < 100) {
      warnings <- c(warnings, i)
    }
  }


  fl <- fsApply(flowset, meanMedianSD) %>% tibble::as_tibble(rownames = "name")

  if(!is.na(channel)) {
    fl <- fl  %>%
    dplyr::select(dplyr::contains(c("name", channel)))}

   name <- fsApply(flowset, function(x) keyword(x)$GUID)
   colnames(name) <- "name"

  if (length(warnings) != 0) {
    warnings <- paste(warnings, collapse = ", ")
    print(paste("Warning: frame(s)", warnings,
                "had less than 100 events in this gate."))
  }

  # Put it all together
   if(exists("atime")){
     flsummary <- dplyr::bind_cols(name, time = time, btime = unname(btime), atime = unname(atime), events = events, conc = conc)
   } else
     flsummary <- dplyr::bind_cols(name, time = time, btime = unname(btime), events = events, conc = conc)
  flsummary <- dplyr::left_join(flsummary, pData(flowset), by = "name")
  # Make rows filename keys
  #rownames(flsummary) <- name
  flsummary <- dplyr::left_join(flsummary, fl, by = "name")
  return(flsummary)
}

#' Generate summary statistics for a flowSet
#'
#' @description Gates a sample to all yeast, then singlet, then doublets.
#' Also calculates singlet to doublet ratio.
#' Returns a list of data frames, e.g. output$singlets, output$doublets, etc.
#'
#'
#' @param flowset the \code{flowSet} to be summarized
#' @param channel \code{character vector} which data channel(s) should be
#' summarized? If excluded (or \code{NA}) all channels will be summarized
#' (default)
#' @param gated \code{boolean} is the data already appropriately gated?
#' @param ploidy \code{character} does the flowSet contain haploid or diploid
#' cells?
#' @param only \code{character} summarize only "singlet", "doublet", or all
#' "yeast" cells, FALSE will return all
#'
#' @return \code{data frame} containing the specified summary statistics of
#' the specified cell populations for each frame
#'
#' @export
#'
#' @examples
#' plate1 <- read.flowSet(path = system.file("extdata", "ss_example",
#' package = "flowTime"), alter.names = TRUE)
#' summarizeFlow(plate1, channel = "FL1.A", gated = TRUE,
#' ploidy = "diploid", only = "yeast")
#'
summarizeFlow <- function(flowset, channel = NA, gated = FALSE,
                          ploidy = FALSE, only = FALSE) {
  # Number of experiments
  n_experiments <- length(flowset)

  # Gate the samples
  if (gated == TRUE) {
    print("Summarizing all events...")
    sum <- flsummary(flowset, channel = channel)
  }
  else {
    if (ploidy == "haploid") {
      if (only == FALSE | only == "yeast") {
        if (exists("yeastGate")){
          print("Gating with haploid yeast gate...")
          subset <- Subset(flowset, yeastGate)
          print("Summarizing all yeast events...")
          sum <- flsummary(subset, channel = channel)
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
          sum <- flsummary(subset, channel = channel)
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
          sum <- flsummary(subset, channel = channel)
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
        print("'only' must be 'singlets','doublets', or 'yeast'")
        stop()
      }
    }
    else if (ploidy == "diploid") {
      if (only == FALSE | only == "yeast") {
        if (exists("yeastGate")){
          print("Gating with diploid yeast gate...")
          subset <- Subset(flowset, yeastGate)
          print("Summarizing all yeast events...")
          sum <- flsummary(subset, channel = channel)
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
          sum <- flsummary(subset, channel = channel)
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
          sum <- flsummary(subset, channel = channel)
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
        print("'only' must be 'singlets','doublets', or 'yeast'")
        stop()
      }
    }
    else {
      print("No ploidy specified. Summarizing all events.")
      sum <- flsummary(flowset, channel = channel)
    }
  }
  return(sum)
}

#' Normalize fluorescence
#'
#' @importFrom plyr "ddply"
#' @description Produces a normalized fluorescence column 'normed'. Expects
#' the 'FL1.A_bs' column to exist or a column to be specified. Has three
#' different methods, version 1 and version 2, described in the script
#'
#' @param frame \code{data frame} of summary statistics to be normalized
#' @param factor_in \code{character vector} containing the varibles to split
#' the data frame by
#' @param method which normalization method to use, 1, 2 or 3.
#' @param column \code{character} the column to apply the normalization to
#'
#' @return \code{data frame} containing the additional normalized variable
#' @details Method 1, the default normalization method, takes the highest
#' point in each dataset grouped by 'factor_in' and normalizes all values in
#' the group by this point. This method is default because it works
#' regardless of whether the data is a time series. Method 2 finds the mean
#' value of all time points with time values less than 0 for each group and
#' normalizes each group by this respective value. Requires a time series
#' with negative time values to work. Method 3 fits a linear model to the
#' pre-zero time points for each groups, infers the y-intercept, and
#' normalizes using this intercept. Method 3 also requires a
#' time series with negative time values to work.
#' @export
#'
#' @examples
#' dat <- read.flowSet(path=system.file("extdata", "tc_example",
#' package = "flowTime"), alter.names = TRUE)
#' annotation <- read.csv(system.file("extdata", "tc_example.csv",
#' package = "flowTime"))
#' adat <- annotateFlowSet(dat, annotation)
#' loadGates(gatesFile = 'C6Gates')
#' dat_sum <- summarizeFlow(adat, ploidy = "diploid", only = "singlets",
#' channel = "FL1.A")
#' dat_sum <- addnorm(dat_sum, c("strain", "treatment"), method = 1,
#' column = "FL1.Amean")
#'
addnorm <- function(frame, factor_in = c("strain", "treatment"),
                    method = 1, column = "FL3.Amean_bs") {
  if ((sum(colnames(frame) == column)) == 0) {
    if ((sum(colnames(frame) == "FL3.A_bs")) == 0) {
      stop("Could not find the background-subtracted values column. \n
           This script requires that there be a column named \n
           FL1.Amean_bs, FL1.A_bs, or the user-defined column using\n
           column='desired-column'")
    } else {
      column <- "FL3.A_bs"
    }
  }
  factors <- which(lapply(frame, class) == "factor")
  frame[,factors] <- sapply(frame[,factors],as.character)
  if (method == 1) {
    # Default normalization method. Takes highest point in dataset grouped
    # by 'factor_in' and sets
    # it to 1, divides all other values by that number. This method is
    # default because it works
    # regardless of whether the data is a time series.
    estimate_0 <- function(x) {
      x[, "normed"] = x[, column]/max(x[, column])
      return(x)
    }
  } else if (method == 2) {
    # Version 2 - takes the mean value of all time points which are less than
    # 0, after grouped by
    # 'factor_in'.  Sets this to the value by which all other data points in
    # that group are divided
    # Therefore, no value is actually '1' except by very rare chance
    # Requires a time series with negative time values to work
    estimate_0 <- function(x) {
      normresult <- x[, column]/mean(x[x$time < 0, column])
      x <- cbind(x, normed = normresult)
      return(x)
    }
  } else if (method == 3) {
    # Version 3 makes a fit line to all pre-zero time points and infers the
    # y-intercept.
    # Requires a time series with negative time values to work.
    estimate_0 <- function(x) {
      prezero_points <- x[x$time < 0, ]
      prezero_fit <- stats::lm(prezero_points[, column] ~ prezero_points[, "time"])
      prezero_intercept <- prezero_fit$coefficients[1]  # intercept
      normresult <- x[, column]/prezero_intercept
      x <- cbind(x, normed = normresult)
      return(x)
    }
  } else {
    stop("You must define version=1, version=2, or version=3)")
  }

  # Check for negative time values
  if (sum(frame$time < 0) == 0) {
    if (method == 2 | method == 3) {
      stop("To use methods 2 or 3, the input data frame must have negative
           time values for each normalized data subset")
    }
  }

  # Run the chosen estimation function and apply it
  frame <- plyr::ddply(frame, factor_in, estimate_0)
  frame[,factors] <- apply(frame[,factors],2,as.factor)
  return(frame)
}


#' Add background subtraction to a summary data frame
#' @description Makes a new column from \code{column} with the background value of
#' a given \code{baseline} control from a chosen identifier column
#' \code{baseline_column} subtracted from the values of \code{column}.
#'
#' @param data the summary data frame of a flowSet (from
#' \code{\link{summarizeFlow}} or
#' \code{\link{flsummary}}) to be used in calculating the background
#' subtracted column
#' @param column the column containing the fluorescent (or other)
#' measurement to be background subtracted
#' @param baseline_column the column containing the identifier of the
#' rows containing background values
#' @param baseline \code{character} the identified or name of representing
#' background fluorescent values
#'
#' @return A summary data frame with an additional column \code{column_bs}
#' containing the background subtracted values
#' @export
#' @importFrom rlang :=
#'
#' @examples
#' dat<-read.flowSet(path=system.file("extdata", "tc_example",
#' package = "flowTime"),alter.names = TRUE)
#' annotation <- read.csv(system.file("extdata", "tc_example.csv",
#' package = "flowTime"))
#' annotation[which(annotation$treatment == 0), 'strain'] <- 'background'
#' adat <- annotateFlowSet(dat, annotation)
#' dat_sum <- summarizeFlow(adat, gated = TRUE,
#' channel = 'FL1.A')
#' dat_sum <- addbs(data = dat_sum, column = FL1.Amean,
#' baseline_column = strain,
#' baseline = "background")
#'
addbs <- function(data, column, baseline_column,
                  baseline = "noYFP") {
      data %>% dplyr::mutate("{{column}}_bs" := {{column}} -
                     mean({{column}}[{{baseline_column}} == baseline]))
}
