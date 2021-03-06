% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summaryAnalysis.R
\name{addnorm}
\alias{addnorm}
\title{Normalize fluorescence}
\usage{
addnorm(
  frame,
  factor_in = c("strain", "treatment"),
  method = 1,
  column = "FL3.Amean_bs"
)
}
\arguments{
\item{frame}{\code{data frame} of summary statistics to be normalized}

\item{factor_in}{\code{character vector} containing the varibles to split
the data frame by}

\item{method}{which normalization method to use, 1, 2 or 3.}

\item{column}{\code{character} the column to apply the normalization to}
}
\value{
\code{data frame} containing the additional normalized variable
}
\description{
Produces a normalized fluorescence column 'normed'. Expects
the 'FL1.A_bs' column to exist or a column to be specified. Has three
different methods, version 1 and version 2, described in the script
}
\details{
Method 1, the default normalization method, takes the highest
point in each dataset grouped by 'factor_in' and normalizes all values in
the group by this point. This method is default because it works
regardless of whether the data is a time series. Method 2 finds the mean
value of all time points with time values less than 0 for each group and
normalizes each group by this respective value. Requires a time series
with negative time values to work. Method 3 fits a linear model to the
pre-zero time points for each groups, infers the y-intercept, and
normalizes using this intercept. Method 3 also requires a
time series with negative time values to work.
}
\examples{
dat <- read.flowSet(path=system.file("extdata", "tc_example",
package = "flowTime"), alter.names = TRUE)
annotation <- read.csv(system.file("extdata", "tc_example.csv",
package = "flowTime"))
adat <- annotateFlowSet(dat, annotation)
loadGates(gatesFile = 'C6Gates')
dat_sum <- summarizeFlow(adat, ploidy = "diploid", only = "singlets",
channel = "FL1.A")
dat_sum <- addnorm(dat_sum, c("strain", "treatment"), method = 1,
column = "FL1.Amean")

}
