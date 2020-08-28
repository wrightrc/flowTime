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
  mdata <- dplyr::left_join(pData(yourFlowSet), annotation_df, by = mergeBy)
  rownames(mdata) <- as.character(mdata[, mergeBy])
  pData(yourFlowSet) <- mdata
  yourFlowSet
}


#' Summary statistic columns for a flow frame
#'
#' @param frame a \code{flowFrame}
#'
#' @return a matrix with a single row for the flow frame and mean, median, and
#' sd columns for each column of the expression measurements
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom tibble "as_tibble"
#'
#' @examples
#' plate1<-read.flowSet(path = system.file("extdata", "ss_example", package =
#' "flowTime"), alter.names = TRUE)
#' meanMedianSD(plate1@frames$A01.fcs)
#' fsApply(plate1, meanMedianSD)
meanMedianSD <- function(frame){
  name <- frame@description$GUID
  frame <- exprs(frame)
  out <- frame %>%
    tibble::as_tibble() %>%
    dplyr::summarise(dplyr::across(
      .fns = list(mean = mean,
                  median = stats::median,
                  sd = stats::sd
                  ),
      .names = "{col}{fn}")) %>% as.matrix
  row.names(out) <- name
  out
}

#' Quality assurance check
#' @description Check whether a flowSet (or a single flowFrame) contains empty
#' values, in which case normalization may fail (divide by zero). This is
#' particularly useful for removing wash wells from a flowSet.
#'
#' @param x \code{flowSet} or \code{flowFrame} to be checked
#' @param threshold \code{flowFrames} with fewer events than this threshold
#' will be identified.
#'
#' @return A vector containing the \code{flowFrames} with fewer events than
#' the threshold.
#' @export
#'
#' @examples
#' plate1<-read.flowSet(path = system.file("extdata", "ss_example", package =
#' "flowTime"), alter.names = TRUE)
#' qaGating(plate1)
#'
qaGating <- function(x, threshold = 100) {
  # Defaults to event count threshold of 100

  print("Running QA...")
  x.class <- class(x)[1]
  if (x.class == "flowFrame") {
    counts <- length(exprs(x[, 1]))
  } else if (x.class == "flowSet") {
    counts <- fsApply(x, length, use.exprs = TRUE)
  } else {
    print("Input must be a flowSet or flowFrame")
  }

  # Find all counts less than 100 (now it's a boolean vector)
  counts.boolean <- counts < threshold
  #positions of those that failed
  counts.failed.position <- grep(TRUE, counts.boolean)

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
#' @param flowframe The \code{flowFrame} for which you would like the initial
#' time
#'
#' @return \code{numeric} time value in minutes
#' @export
#'
#' @examples
#' plate1<-read.flowSet(path = system.file("extdata", "ss_example", package =
#' "flowTime"),alter.names = TRUE)
#' getTime(plate1$A01.fcs)
getTime <- function(flowframe) {
  time_raw <- as.numeric(unlist(strsplit(keyword(flowframe)$`$BTIM`,
                                         split = ":")))
  time <- time_raw[1] * 60 + time_raw[2] + time_raw[3]/60 + time_raw[4]/6000
  return(time)
}


#' Read FCS files from set of plates
#'
#' @description Reads all folders within the specified path containing the
#' specified pattern in the folder names. Each folder contains a set a plate
#' of FCS files. These folders typically make up a whole experiment. Plates
#' are numbered according to the standard lexicographical ordering of your
#' operating system.
#'
#' @param path The path to search for folders containing FCS files
#' @param pattern The \link[base]{regex} pattern used to identify
#' the folders of FCS files to be read
#' @param ... Additional arguments passed to read.flowSet.
#' Note that `alter.names` is forced to be TRUE in this implementation.
#'
#' @return A single flowSet containing all FCS files within the
#' identified folders. The index of each folder in the list according
#' to lexicographical ordering (1,2,...) is prepended to the sampleNames.
#'
#' @export
#'
#' @examples
#' # Read in both of the example data sets as a single flowSet
#' plate1<-read.plateSet(path = system.file("extdata", package = "flowTime"),
#' pattern = "")
#'
read.plateSet <- function(path = getwd(), pattern = "", ...){
  files <- grep(pattern = pattern, x = list.dirs(path, recursive = FALSE,
                                                 full.names = FALSE),
                value = TRUE)
  for(file in files) {
    plate <- read.flowSet(path = paste(path, file, sep = "/"), alter.names = TRUE, ...)
    plate_num <- which(files == file)
    sampleNames(plate) <- paste0(plate_num, sampleNames(plate))
    pData(plate)$name <- sampleNames(plate)
    pData(plate)$folder <- file
    if('flow_set' %in% ls()) flow_set <- rbind2(flow_set, plate) else flow_set <- plate
  }
  return(flow_set)
}
