## ---- echo = FALSE, message=FALSE---------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", tidy = TRUE)
library(flowCore)
library(flowTime)
library(ggplot2)

## -----------------------------------------------------------------------------
plate1<-read.flowSet(path=system.file("extdata", "ss_example", 
                        package = "flowTime"), alter.names = TRUE)
# add plate numbers to the sampleNames
sampleNames(plate1)<-paste("1_", sampleNames(plate1), sep = "")
dat<-plate1

## ---- eval = F----------------------------------------------------------------
#  plate2 <- read.flowSet(path = paste(experiment, "_2/", sep = ""),
#                          alter.names = TRUE)
#  sampleNames(plate2) <- paste("2_", sampleNames(plate2), sep = "")
#  dat <- rbind2(plate1, plate2)

## -----------------------------------------------------------------------------
annotation <- read.csv(system.file("extdata", "ss_example.csv", 
                        package = "flowTime"))
head(annotation)
sampleNames(dat) 
sampleNames(dat) == annotation$name 

## ---- eval = F----------------------------------------------------------------
#  annotation <- cbind(annotation, 'name' =  sampleNames(dat))
#  # or
#  annotation <- createAnnotation(yourFlowSet = dat)
#  write.csv(annotation, file = 'path/to/yourAnnotation.csv')

## -----------------------------------------------------------------------------
adat <- annotateFlowSet(yourFlowSet = dat, annotation_df = annotation, 
                        mergeBy = 'name')
head(rownames(pData(adat)))
head(pData(adat))

## ---- eval = F----------------------------------------------------------------
#  write.flowSet(adat, outdir = 'your/favorite/directory')
#  
#  # Read the flowSet with the saved experimental meta data
#  read.flowSet('flowSet folder', path = 'your/flow/directory',
#                          phenoData = 'annotation.txt', alter.names = TRUE)

## ---- fig.width = 4, fig.height = 4-------------------------------------------
loadGates() # use the default included gateSet
dat.SS <- steadyState(flowset = adat, ploidy = 'diploid', only = 'singlets')

p <- ggplot(dat.SS, aes(x = as.factor(treatment), y = FL2.A, fill = AFB)) + 
  geom_boxplot(outlier.shape = NA) + facet_grid(IAA~AFB) + 
  theme_classic(base_family = 'Arial', base_size = 16) + ylim(c(-1000,10000)) +
  xlab(expression(paste('Auxin (',mu,'M)',sep = ""))) + 
  ylab('Fluorescence (AU)') + theme(legend.position="none")
p

