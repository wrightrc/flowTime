## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", tidy = TRUE)
library(flowCore)
library(flowTime)
library(ggplot2)

## -----------------------------------------------------------------------------
plate1<-read.flowSet(path=system.file("extdata", "tc_example", package = "flowTime"), alter.names = TRUE)
# add plate numbers to the sampleNames, in this example we have already done this
# step
# sampleNames(plate1)<-paste("1_",sampleNames(plate1),sep="")
dat<-plate1

## ---- eval = F----------------------------------------------------------------
#  plate2<-read.flowSet(path = paste(experiment,"_2/",sep=""), alter.names = TRUE)
#  sampleNames(plate2)<-paste("2_", sampleNames(plate2), sep = "")
#  dat<-rbind2(plate1, plate2)

## -----------------------------------------------------------------------------
annotation <- read.csv(system.file("extdata", "tc_example.csv", package = 
                                     "flowTime"))

## ---- eval = F----------------------------------------------------------------
#  sampleNames(dat) # view the sample names
#  sampleNames(dat) == annotation$id
#  # Replace 'id' with the unique identifier column to test,
#  # if this column is identical to the sample names of your flowset.
#  annotation <- cbind(annotation, 'names' =  sampleNames(dat))
#  # If the sampleNames and unique identifiers are in the correct order
#  # this command will add the sampleNames as the identifier.

## -----------------------------------------------------------------------------
adat <- annotateFlowSet(dat, annotation)
head(rownames(pData(adat)))
head(pData(adat))

## ---- eval = F----------------------------------------------------------------
#  write.flowSet(adat, outdir = 'your/favorite/directory')
#  read.flowSet('flowSet folder', path = 'your/flow/directory',
#               phenoData = 'annotation.txt', alter.names = TRUE)

## ---- fig.width= 7------------------------------------------------------------
#load the gate set for BD Accuri C6 cytometer
loadGates(gatesFile = 'C6Gates.RData', path = system.file("extdata", package = "flowTime"))
dat_sum <- summarizeFlow(adat, ploidy = 'diploid', 
                         only = 'singlets',channel = 'FL1.A')

qplot(x = time, y= FL1.Amean, data = dat_sum, linetype = factor(treatment)) + 
  geom_line() + xlab('Time post Auxin addition (min)') + 
  ylab('Reporter Fluorescence (AU)') + 
  scale_color_discrete(name=expression(paste("Auxin (",mu,"M)",sep = ""))) + 
  theme_classic(base_size = 14, base_family = 'Arial')

