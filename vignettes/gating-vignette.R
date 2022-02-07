## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", tidy = TRUE)
library(flowCore)
library(flowTime)
library(ggplot2)
library(dplyr)

## -----------------------------------------------------------------------------
data <- read.flowSet(path=system.file("extdata", "ss_example", 
                        package = "flowTime"), alter.names = TRUE)
annotation <- read.csv(system.file("extdata", "ss_example.csv", 
                        package = "flowTime"))
adat <- annotateFlowSet(yourFlowSet = data, annotation_df = annotation, 
                        mergeBy = 'name')

## -----------------------------------------------------------------------------
#library(BiocManager)
#BiocManager::install("openCyto")
#BiocManager::install("ggcyto")

library("openCyto")
library("ggcyto")
library("flowClust")

#vignette("flowWorkspace-Introduction", "flowWorkspace")
#vignette('HowToAutoGating', package = "openCyto")


## -----------------------------------------------------------------------------
autoplot(data[21:24], x = "FSC-A", y = "SSC-A") 

## -----------------------------------------------------------------------------
Debris <- gate_quantile(fr = data[[3]], channel = "FSC.A", probs = 0.99, filterId = "Debris")
autoplot(data[[3]], x = "FSC-A") + geom_gate(Debris)

## -----------------------------------------------------------------------------
autoplot(data[21:24], x = "FSC-A", "SSC-A") + geom_gate(Debris)
toTable(summary(filter(data[c(3,21:24)], !Debris)))

## -----------------------------------------------------------------------------
#Initialize the single frame
data.1frame <- data[[1]]
#fill the single frame with the exprs data from each frame 
# in the flow set
exprs(data.1frame) <- fsApply(data, function(x) {
  x <- exprs(x)
  return(x)
})

autoplot(data.1frame, x = "FSC-A", "SSC-A")

## -----------------------------------------------------------------------------
autoplot(data.1frame, x = "FSC-A", "SSC-A") + scale_x_logicle() + scale_y_logicle()

## -----------------------------------------------------------------------------
chnls <- c("FSC.A", "SSC.A", "FSC.H", "SSC.H")
trans <- estimateLogicle(data.1frame, channels = chnls)
inv.trans <- inverseLogicleTransform(trans)
data.1frame <- transform(data.1frame, trans)
autoplot(data.1frame, x = "FSC-A", "SSC-A")

## -----------------------------------------------------------------------------
yeast <- gate_flowClust_2d(data.1frame, xChannel = "FSC.A", 
                           yChannel =  "SSC.A", K = 1, 
                           quantile = .95, min = c(0,0))
autoplot(data.1frame, x = "FSC-A", y = "SSC-A") + geom_gate(yeast)


## -----------------------------------------------------------------------------
yeast <- transform(yeast, inv.trans)
data.1frame <- transform(data.1frame, inv.trans)

## -----------------------------------------------------------------------------
autoplot(data[c(1, 8, 16, 24, 32)], "FSC.A","SSC.A") + 
  geom_gate(yeast)
#invisible(capture.output( 
  # we have to use this to prevent summary from printing
  f<- summary(filter(data, yeast))#))
# Now we can print our summary as a table
toTable(f)

## -----------------------------------------------------------------------------
autoplot(Subset(data.1frame, yeast), "FSC-A", "FSC-H")
library(flowStats)
chnl <- c("FSC-A", "FSC-H")
singlets <- gate_singlet(x = Subset(data.1frame, yeast), area = "FSC.A",
                         height = "FSC.H", prediction_level = 0.999, maxit = 20)
autoplot(Subset(data.1frame, yeast), "FSC-A", "FSC-H") + geom_gate(singlets)

## -----------------------------------------------------------------------------
autoplot(data[c(1:4, 29:32)], x = "FSC-A", y = "FSC-H") + 
  geom_gate(singlets) + facet_wrap("name", ncol = 4)

autoplot(Subset(data[c(1:4, 29:32)], yeast & singlets), x = "FL1-A") + 
  facet_wrap("name", ncol = 4) 

## -----------------------------------------------------------------------------
invisible(capture.output(
  d <- summary(filter(data, yeast & singlets))))
(e <- toTable(d))
e <- left_join(e, pData(data), by = c("sample" = "name"))
ggplot(data = e, mapping = aes(x = as.factor(sample), y = percent)) + geom_point()

## ---- eval = FALSE------------------------------------------------------------
#  data <- Subset(data, yeast & singlets)
#  data_sum <- summarizeFlow(data, channel = c("FL1.A", "FL4.A"), gated = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  saveGates(yeastGate = yeast, dipsingletGate = singlets, fileName = "PSB_Accuri_W303.RData")
#  loadGates(gatesFile = "PSB_Accuri_W303.RData")
#  data_sum <- summarizeFlow(data, channel = c("FL1.A", "FL4.A"), ploidy = "diploid", only = "singlets")

