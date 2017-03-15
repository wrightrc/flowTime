
#' Get summary statistics for fluorescence or other data channels of a flowSet
#'
#' @import plyr
#' @import flowCore
#'
#' @param flowset the \code{flowSet} to create summary statistics for
#' @param channel \code{character} the data channel to summarize
#' @param moments \code{boolean} if \code{TRUE} then split each frame into
#' early, middle, and late events
#'
#' @return A \code{data frame} containing summary statistics (mean, median,
#' SD) for the specified fluorescent channel and time moments of the flowSet.
#' @export
#'
#' @examples plate1 <- read.flowSet(path = system.file("extdata",
#' "ss_example", package = "flowTime"), alter.names = TRUE)
#' flsummary(plate1)
flsummary <- function(flowset, channel = "FL3.A", moments = FALSE) {
  # Number of cells (experiments) in the flowSet
  n_experiment <- length(flowset)

  # Initialize empty matrices/data frames to increase efficiency
  warnings <- c()

  if (moments == TRUE) {
    requireNamespace(moments)
  }

  # Get time of each frame in minutes of the day
  btime_raw <- fsApply(flowset, function(x)
    as.numeric(unlist(strsplit(keyword(x)$`$BTIM`, split = ":"))))
  btime <- apply(btime_raw, 1, function(x) x[1] * 60 + x[2] + x[3]/60 +
                   x[4]/6000)
  time <- btime - min(btime)

  # Acquisition time - how long it took to take the sample, in seconds
  atime <- fsApply(flowset, function(x)
    as.numeric(keyword(x)$`#ACQUISITIONTIMEMILLI`)/1000)

  events <- fsApply(flowset, function(x) length(x[, 1]), use.exprs = TRUE)
  uL <- fsApply(flowset, function(x) as.integer(keyword(x)$`$VOL`)/1000)
  conc <- events/uL

  for (i in 1:n_experiment) {
    if (events[i] < 100) {
      warnings <- c(warnings, i)
    }
  }

  fl_mean <- fsApply(flowset, function(x) base::mean(x[, channel]),
                     use.exprs = TRUE)
  fl_median <- fsApply(flowset, function(x) stats::median(x[, channel]),
                       use.exprs = TRUE)
  fl_sd <- fsApply(flowset, function(x) stats::sd(x[, channel]), use.exprs = TRUE)
  fl <- data.frame(fl_mean, fl_median, fl_sd)
  colnames(fl) <- paste(channel, c("mean", "median", "sd"), sep = "")

  # Do we want the first few moments?
  if (moments == TRUE) {
    requireNamespace(moments)
    fl_var <- data.frame(fsApply(flowset, function(x)
      stats::var(x[, channel]), use.exprs = TRUE))
    fl_skew <- data.frame(fsApply(flowset, function(x)
      moments::skewness(x[, channel]), use.exprs = TRUE))
    fl_kurt <- data.frame(fsApply(flowset, function(x)
      moments::kurtosis(x[, channel]), use.exprs = TRUE))
    fl_moments <- data.frame(fl_var, fl_skew, fl_kurt)
    colnames(fl_moments) <- paste(channel, c("var", "skew", "kurt"), sep = "")
    fl <- cbind(fl, fl_moments)
  }

  name <- fsApply(flowset, function(x) keyword(x)$GUID)
  colnames(name) <- "name"

  if (length(warnings) != 0) {
    warnings <- paste(warnings, collapse = ", ")
    print(paste("Warning: frame(s)", warnings,
                "had less than 100 events in this gate."))
  }

  # Put it all together
  flsummary <- cbind(time, btime, atime, events, conc, fl, name)

  # Make rows filename keys
  rownames(flsummary) <- name
  flsummary <- join(flsummary, pData(flowset), by = "name")
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
#' @param channel \code{character} which data channel should be summarized
#' @param gated \code{boolean} is the data already appropriately gated?
#' @param ploidy \code{character} does the flowSet contain haploid or diploid
#' cells?
#' @param moments \code{boolean} split the data into early, middle, and late
#' moments?
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
#' ploidy = "diploid", moments = FALSE, only = "yeast")
#'
summarizeFlow <- function(flowset, channel = "FL1.A", gated = FALSE,
                          ploidy = FALSE, moments = FALSE, only = FALSE) {

  if(!exists(c("yeastGate", "hapsingletGate", "hapdoubletGate",
               "dipsingletGate", "dipdoubletGate"))) loadGates()
  yeastGate <- get("yeastGate")
  hapsingletGate <- get("hapsingletGate")
  hapdoubletGate <- get("hapdoubletGate")
  dipsingletGate <- get("dipsingletGate")
  dipdoubletGate <- get("dipdoubletGate")

  # Number of experiments
  n_experiments <- length(flowset)


  # Gate the samples
  if (gated == TRUE) {
    print("Summarizing all events...")
    flowsum <- flsummary(flowset, channel = channel, moments = moments)
    return(flowsum)
  }
  else {
    if (!exists(c("yeastGate", "hapsingletGate", "hapdoubletGate",
                  "dipsingletGate", "dipdoubletGate")))
      loadGates()
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

    if (only == FALSE) {
      # Normalize and summarize each subset
      print("Summarizing all yeast events...")
      yeastsum <- flsummary(yeast, channel = channel, moments = moments)
      print("Summarizing doublets events...")
      doubletsum <- flsummary(doublets, channel = channel,
                              moments = moments)

      print("Summarizing singlets events...")
      singletsum <- flsummary(singlets, channel = channel,
                              moments = moments)
    }
    else {
      if (only == "singlets") {
        print("Summarizing singlets events...")
        singletsum <- flsummary(singlets, channel = channel,
                                moments = moments)
        return(singletsum)
      }
        else if (only == "doublets") {
        print("Summarizing doublets events...")
        doubletsum <- flsummary(doublets, channel = channel,
                                moments = moments)
        return(doubletsum)
      }
        else if (only == "yeast") {
        print("Summarizing all yeast events...")
        yeastsum <- flsummary(yeast, channel = channel,
                              moments = moments)
        return(yeastsum)
      }
      else {
        print("'only' must be 'singlets','doublets', or 'yeast'")
        stop()
      }
    }
  }
  summary_list <- list(yeast = yeastsum, singlets = singletsum,
                       doublets = doubletsum)
  return(summary_list)
}

#' Normalize fluorescence
#' @description Produces a normalized fluorescence column 'normed'. Expects
#' the 'FL1.A_bs' column to exist or a column to be specified. Has two
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
#' with negative time values to work. Version 3 fits a linear model to the
#' pre-zero time points for each groups,  infers the y-intercept, and
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
  frame <- ddply(frame, factor_in, estimate_0)
  frame[,factors] <- apply(frame[,factors],2,as.factor)
  return(frame)
}


#' Add background subtraction to a summary data frame
#' @description Subtracts the background fluorescence of a given control
#' strain from the chosen column.
#'
#' @param flowData the summary data frame of flowSet to be background
#' subtracted
#' @param column the column containing the fluorescent measurement to be
#' background subtracted
#' @param baseline_column the column containing the name of the strain
#' representing background fluorescent values
#' @param baseline \code{character} the name of the strain representing
#' background fluorescent values
#'
#' @return A summary data frame with an additional column "column_bs"
#' containing the background subtracted fluorescent values
#' @export
#'
#' @examples
#' dat<-read.flowSet(path=system.file("extdata", "tc_example",
#' package = "flowTime"),alter.names = TRUE)
#' annotation <- read.csv(system.file("extdata", "tc_example.csv",
#' package = "flowTime"))
#' annotation[which(annotation$treatment == 0), 'strain'] <- 'background'
#' adat <- annotateFlowSet(dat, annotation)
#' loadGates(gatesFile = 'C6Gates')
#' dat_sum <- summarizeFlow(adat, ploidy = 'diploid', only = 'singlets',
#' channel = 'FL1.A')
#' dat_sum <- addbs(dat_sum, column = "FL1.Amean", baseline = "background")
#'
addbs <- function(flowData, column = "FL3.Amean", baseline_column = "strain", baseline = "noYFP") {
  flowData[, paste(column, "_bs", sep = "")] <- flowData[, column] -
    mean(subset(flowData, baseline_column == baseline)[,column])
  return(flowData)
}
