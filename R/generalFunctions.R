#' Quality assurance check
#' @description Check whether a flowSet (or a single flowFrame) contains empty values, in which case normalization may fail (divide by zero). This is particularly useful for removing wash wells from a flowSet.
#'
#' @param x \code{flowSet} or \code{flowFrame} to be checked
#' @param threshold \code{flowFrames} with fewer events than this threshold will be identified.
#'
#' @return A vector containing the \code{flowFrames} with fewer events than the threshold.
#' @export
#'
#' @examples
#' plate1<-read.flowSet(path = system.file("extdata", "ss_example", package = "flowTime"), alter.names = TRUE)
#' qa.gating(plate1)
#'
qa.gating <- function(x, threshold = 100) {
  # Defaults to event count threshold of 100

  print("Running QA...")
  x.class <- class(x)[1]
  if (x.class == "flowFrame") {
    counts <- length(exprs(x[, 1]))
  } else if (x.class == "flowSet") {
    counts <- fsApply(x, length, use.exprs = T)
  } else {
    print("Input must be a flowSet or flowFrame")
  }

  # Find all counts less than 100 (now it's a boolean vector)
  counts.boolean <- counts < threshold
  counts.failed.position <- grep(TRUE, counts.boolean)  #positions of those that failed

  # Did we get any failed counts?  If so, return the position
  counts.failed <- length(counts.failed.position) != 0
  if (counts.failed) {
    print("QA resulted in 1 or more warnings.")
    return(counts.failed.position)
  } else {
    print("QA succeeded")
    return(FALSE)
  }
}

#' Get the time at which at flowFrame began collection
#'
#' @param flowframe The \code{flowFrame} for which you would like the initial time
#'
#' @return \code{numeric} time value in minutes
#' @export
#'
#' @examples
#' plate1<-read.flowSet(path = system.file("extdata", "ss_example", package = "flowTime"),alter.names = TRUE)
#' get_time(plate1$A01.fcs)
get_time <- function(flowframe) {
  time_raw <- as.numeric(unlist(strsplit(keyword(flowframe)$`$BTIM`, split = ":")))
  time <- time_raw[1] * 60 + time_raw[2] + time_raw[3]/60 + time_raw[4]/6000
  return(time)
}
