### fl1transform:
### Normalizes FL1.A values in a flowset/flowframe to FSC.A values
### Should control for (at least some) non-linearity in the values
### Also makes FL1.A proportional to fluorescence/cell volume
### Used to multiply by the mean FSC.A value, but this is probably
### statistically questionable. Now multiplies by a constant (10000)
### simply to keep the values in the integer range.
###
### If you specify transform="log", it will simply do a log transform
### to FL1.A instead.

#' Title
#'
#' @param x
#' @param transform
#'
#' @return
#' @export
#'
#' @examples
fl1transform <- function(x,transform=F) {
  # Default scaling is 10^4 and is solely for making results human-readable

  # Handle both flowFrames and flowSets
  x.class <- class(x)[1]
  if (x.class=="flowFrame") {
    #stop("This is a flowFrame, not a flowSet")
    return(transform(x,FL1.A=FL1.A/FSC.A*10^4))
  }

  # Protect the input from the modifications
  x <- x[seq(along=x)]

  # Do QA.  Reject all frames with no cells
  qa.result <- qa.gating(x,threshold=1)

  # Remove all 0-valued fluorescence results.
  # These are very likely to be artifacts and screw up some transforms
  print("Removing 0-valued fluorescence outliers")
  x <- Subset(x,rectangleGate(FL1.A=c(0.001,Inf)))

  # Transformation setup
  trans <- function(x) {
    if (transform == "fscanorm") {
      x <- transform(x,FL1.A=FL1.A/FSC.A*10^4)
    } else if (transform == "log") {
      x <- transform(x,FL1.A=log(FL1.A))
    } else if (transform == F) {
      x <- x # do nothing.  Is this necessary?
    } else {
      stop("No legitimate transform set.  Use transform=\"log\" or transform=\"fscanorm\".")
    }
  }

  if (!qa.result) {
    x <- trans(x)
    #		x <- transform(x,FL1.A=FL1.A/FSC.A*10^4)
  } else {
    # For loop is inefficient, switch this out for fsApply while maintaining all cells
    for (i in qa.result) {
      x[[i]] <- trans(x)
      #			x[[i]] <- transform(x[[i]],FL1.A=FL1.A/FSC.A*10^4)
    }
    #		x <- fsApply(x,transform,FL1.A=FL1.A/FSC.A*10^4)
    cat(paste(
      "### Too few cells at this gating level for frame(s) \n### ",
      paste(qa.result,collapse=", "),
      ".\n### These frames were not normalized.\n\n",sep=""))
  }

  return(x)
}

### flsummary:
### Get summary statistics for fluorescence, other data

#' Title
#'
#' @param flowset
#' @param channel
#' @param moments
#' @param split
#' @param transform
#'
#' @return
#' @export
#'
#' @examples
flsummary <- function(flowset,channel="FL3.A",moments=F,split=F,transform=F) {
  # Number of cells (experiments) in the flowSet
  n_experiment <- length(flowset)

  # Initialize empty matrices/data frames to increase efficiency
  warnings <- c()

  if (moments == T) {
    library(moments)
  }

  # Get time of each frame in minutes of the day
  btime_raw <- fsApply(flowset,function(x)as.numeric(unlist(strsplit(keyword(x)$`$BTIM`,split=":"))))
  btime <-  apply(btime_raw,1,function(x)x[1]*60+x[2]+x[3]/60+x[4]/6000)
  time <- btime-min(btime)

  # Acquisition time - how long it took to take the sample, in seconds
  atime <- fsApply(flowset,function(x)as.numeric(keyword(x)$`#ACQUISITIONTIMEMILLI`)/1000)

  events <- fsApply(flowset,function(x)length(x[,1]),use.exprs=T)
  uL <- fsApply(flowset,function(x)as.integer(keyword(x)$`$VOL`)/1000)
  conc <- events/uL

  for (i in 1:n_experiment) {
    if (events[i] < 100) {
      warnings <- c(warnings,i)
    }
  }

  fl_mean <- fsApply(flowset,function(x)mean(x[,channel]),use.exprs=T)
  fl_median <- fsApply(flowset,function(x)median(x[,channel]),use.exprs=T)
  fl_sd <- fsApply(flowset,function(x)sd(x[,channel]),use.exprs=T)
  fl <- data.frame(fl_mean,fl_median,fl_sd)
  colnames(fl) <- paste(channel,c("mean","median","sd"),sep="")

  # Do we want mean fl values for data split into 4 evenly sized chunks?
  if (split==T) {
    split_table <- fsApply(flowset,splitFrame)
    split_table <- data.frame(matrix(unlist(split_table),ncol=4,byrow=T))
    colnames(split_table) <- paste("split",1:4,sep="")
    fl <- cbind(fl,split_table)
  }

  # Do we want the first few moments?
  if (moments == T) {
    require(moments)
    fl_var <- data.frame(fsApply(flowset,function(x)var(x[,channel]),use.exprs=T))
    fl_skew <- data.frame(fsApply(flowset,function(x)skewness(x[,channel]),use.exprs=T))
    fl_kurt <- data.frame(fsApply(flowset,function(x)kurtosis(x[,channel]),use.exprs=T))
    fl_moments <- data.frame(fl_var,fl_skew,fl_kurt)
    colnames(fl_moments) <- paste(channel,c("var","skew","kurt"),sep="")
    fl <- cbind(fl,fl_moments)
  }

  name <- fsApply(flowset,function(x) keyword(x)$GUID)
  colnames(name) <- "name"

  if (length(warnings) != 0) {
    warnings <- paste(warnings,collapse=", ")
    print(paste("Warning: frame(s)",warnings,"had less than 100 events in this gate."))
  }

  # Put it all together
  flsummary <- cbind(time,btime,atime,events,conc,fl,name)

  # Make rows filename keys
  rownames(flsummary) <- name
  flsummary <- join(flsummary,pData(flowset), by = 'name')

  # Rename the 'mean', 'median', and 'sd' columns to reflect transformations done or channel used.
  # 'FL1.A' = no transformation, 'FL1_FSC' = "fsacanorm", 'log' = "log"
  flsummary <- renameflcols(flsummary,channel=channel,transform=transform)

  return(flsummary)
}

#' rename flsummary data
#' keeping it separate for ease of use
#' probably slow due to for loop
#'
#' @param x
#' @param channel
#' @param transform
#'
#' @return
#' @export
#'
#' @examples
renameflcols <- function(x,channel="FL1.A",transform=F) {
  cols <- c("mean","median","sd")
  if (transform!=F) {
    if (transform=="fscanorm") {
      tname <- "FL1_FSC"
    } else if (transform=="log") {
      tname <- "log"
    } else {
      stop("invalid transform")
    }
  } else {
    return(x)
  }
  for (i in cols) {
    colnames(x)[grep(i,paste(channel,colnames(x),sep=""))] <- paste(tname,i,sep="")
  }
  return(x)
}

#' summary.cyt:
#' Gates a sample to all yeast, then singlet, then doublets
#' Does the work of singletsummary.cyt,doubletsummary.cyt,yeastsummary.cyt
#' Also calculates singlet to doublet ratio
#' Returns a list of data frames, e.g. output$singlets, output$doublets, etc.
#'
#'
#' @param flowset
#' @param transform
#' @param channel
#' @param gated is the data already appropriately gated?
#' @param ploidy
#' @param moments
#' @param split
#' @param only
#'
#' @return
#' @export
#'
#' @examples
summary.cyt <- function(
  flowset,
  transform=F,
  channel="FL1.A",
  gated = F,
  ploidy=F,
  moments=F,
  split=F,
  only=F) {

  # Number of experiments
  n_experiments <- length(flowset)

  # If using channel="FSC.A", don't use fscanorm
  if (channel=="FSC.A"&transform=="fscanorm") {
    print("Channel FSC.A selected with no transform= setting set.")
    print("Defaulting to no transform (set transform=\"log\" for log transform)")
    transform=F
  }

  # Transform FL1.A
  if (transform != F) {
    print(paste("Transforming FL1.A using",
                transform,
                "transform..."
    )
    )
    flowset <- fl1transform(flowset,transform=transform)

  }


  # Gate the samples
  ### Gate the samples
  if (gated == F) {
    if (!exists(c('yeastGate','hapsingletGate','hapdoubletGate','dipsingletGate','dipdoubletGate'))) loadGates()
    if (ploidy=="haploid") {
      print("Gating with haploid gates...")
      yeast <- Subset(flowset,yeastGate)
      singlets <- Subset(yeast,hapsingletGate)
      doublets <- Subset(yeast,hapdoubletGate)
    } else if (ploidy=="diploid") {
      print("Gating with diploid gates...")
      yeast <- Subset(flowset,yeastGate)
      singlets <- Subset(yeast,dipsingletGate)
      doublets <- Subset(yeast,dipdoubletGate)
    } else {
      stop('Error: You must define ploidy="haploid" or ploidy="diploid"')
    }
  }
  if (only==F) {
    # Normalize and summarize each subset
    print("Summarizing all yeast events...")
    yeastsum <- flsummary(yeast,channel=channel,moments=moments,split=split,transform=transform)

    print("Summarizing doublets events...")
    doubletsum <- flsummary(doublets,channel=channel,moments=moments,split=split,transform=transform)

    print("Summarizing singlets events...")
    singletsum <- flsummary(singlets,channel=channel,moments=moments,split=split,transform=transform)
  } else {
    if (only=="singlets") {
      print("Summarizing singlets events...")
      singletsum <- flsummary(singlets,channel=channel,moments=moments,split=split,transform=transform)
      return(singletsum)
    } else if (only=="doublets") {
      print("Summarizing doublets events...")
      doubletsum <- flsummary(doublets,channel=channel,moments=moments,split=split,transform=transform)
      return(doubletsum)
    } else if (only=="yeast") {
      print("Summarizing all yeast events...")
      yeastsum <- flsummary(yeast,channel=channel,moments=moments,split=split,transform=transform)
      return(yeastsum)
    } else {
      print("'only' must be 'singlets','doublets', or 'yeast'")
      stop()
    }
  }

  summary_list <- list(yeast=yeastsum,singlets=singletsum,doublets=doubletsum)
  return(summary_list)
}

# Produces a normalized fluorescence column 'normed'
# Expects the 'FL1.A_bs' column to exist (not hard to extend to others/make it user selectable)
# Has two different methods, version 1 and version 2, described in the script
#' Title
#'
#' @param frame
#' @param factor_in
#' @param method
#' @param column
#'
#' @return
#' @export
#'
#' @examples
addnorm <- function(frame,factor_in=c("strain","treatment"),method=1,column="FL3.Amean_bs") {
  library(plyr)
  if ( (sum(colnames(frame)==column)) == 0 ) {
    if( (sum(colnames(frame)=="FL3.A_bs")) == 0 ) {
      stop("Could not find the background-subtracted values column. \
           This script requires that there be a column named \
           FL1.Amean_bs, FL1.A_bs, or the user-defined column using\
           column='desired-column'")
    } else {
      column <- "FL3.A_bs"
    }
  }

  if (method==1) {
    # Default normalization method. Takes highest point in dataset grouped by 'factor_in' and sets it to 1,
    # divides all other values by that number. This method is default because it works regardless of
    # whether the data is a time series.
    estimate_0 <- function(x) {
      x[,"normed"]=x[,column]/max(x[,column])
      return(x)
    }
  } else if (method == 2) {
    # Version 2 - takes the mean value of all time points which are less than 0, after grouped by 'factor_in'.
    # Sets this to the value by which all other data points in that group are divided
    # Therefore, no value is actually '1' except by very rare chance
    # Requires a time series with negative time values to work
    estimate_0 <- function(x) {
      normresult <- x[,column]/mean(x[x$time<0,column])
      x <- cbind(x,normed=normresult)
      return(x)
    }
  } else if (method == 3) {
    # Version 3 makes a fit line to all pre-zero time points and infers the y-intercept
    # Requires a time series with negative time values to work
    estimate_0 <- function(x) {
      prezero_points <- x[x$time<0,]
      prezero_fit <- lm(prezero_points[,column]~prezero_points[,"time"])
      prezero_intercept <- prezero_fit$coefficients[1] # intercept
      normresult <- x[,column]/prezero_intercept
      x <- cbind(x,normed=normresult)
      return(x)
    }
  } else {
    stop("You must define version=1, version=2, or version=3)")
  }

  # Check for negative time values
  if (sum(frame$time<0)==0) {
    if (method==2|method==3) {
      stop("To use methods 2 or 3, the input data frame must have negative time values for each normalized data subset")
    }
  }

  # Run the chosen estimation function and apply it
  frame <- ddply(frame,factor_in,estimate_0)
  return(frame)
}


addbs <- function(frame,column="FL3.Amean",baseline="noYFP") {
  frame[,paste(column,"_bs",sep="")] <- frame[,column]-mean(subset(frame,strain==baseline)[,column])
  return(frame)
}

