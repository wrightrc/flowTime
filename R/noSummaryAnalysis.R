###########################################################
# R-functions for the analysis of Flow Cytometry data
# using full datasets for each point without prior summary
###########################################################

#' Create an annotation dataframe
#'
#' @param yourFlowSet
#'
#' @return annotation_df
#' @export
#'
#' @examples
createAnnotation <- function(yourFlowSet){
  annotation_df <- data.frame(name = sampleNames(yourFlowSet))
}

# Annotate the flowSet  --------
#' @description  Add annotations to a flowSets phenoData and plate numbers, strain names, and treatment also set T0
#' @param yourFlowSet, a flowSet with sampleNames of the format 'plate#_Well', we typically use the following code chunk to read data from individual plates as exported from BD Accuri C6 software.
#' plate0<-read.flowSet(path=paste(experiment,"/",sep=""),alter.names=T)
#' sampleNames(plate0)<-paste("0_",sampleNames(plate0),sep="")
#' plate1<-read.flowSet(path=paste(experiment,"_1/",sep=""),alter.names=T)
#' sampleNames(plate1)<-paste("1_",sampleNames(plate1),sep="")
#' fullSet<-rbind2(plate0,plate1)
#' @param annotation_df, A data frame with columns "well", "strain", "treatment", containing all of the wells in the flowset labeled with the strain and treatment in that well.
#' Get a vector of all of the wells in your flowset with:
#' well <- fsApply(testSet,function(x) strsplit(keyword(x)$GUID,".fcs")[[1]])
#' Assemble corresponding strain and treatment vectors and then cbind them together.
#' @return An annotated flowSet
#'
#' @examples
annotateFlowSet <- function(yourFlowSet, annotation_df, mergeBy = 'name'){
  pData(yourFlowSet)$name <- sampleNames(yourFlowSet)
  mdata <- join(pData(yourFlowSet),annotation_df, by = mergeBy)
  rownames(mdata)<-as.character(mdata[,mergeBy])
  pData(yourFlowSet) <- mdata
  yourFlowSet
}


# Set T0 ------------------------------------------------------------------
#
# #find time of the beginning of acquisition for each well in minutes
# btime_raw <- fsApply(yourFlowSet,function(x)as.numeric(unlist(strsplit(keyword(x)$`$BTIM`,split=":"))))
# btime <-  apply(btime_raw,1,function(x)x[1]*60+x[2]+x[3]/60+x[4]/6000)
# #Set the cell-Time relative to well-Time
# #set time relative to T0
# time <- btime-btime[T0]
# yourFlowSet<-fsApply(yourFlowSet, function(well) {
#   minTime<-min(exprs(well)[,'Time'])
#   timeStep<-as.numeric(keyword(well)$`$TIMESTEP`)
#   wellTime<-pData(yourFlowSet)[identifier(well),'time']
#   exprs(well)[,'Time']<-apply(exprs(well$Time),1,function(cellTime) {
#     (cellTime-minTime)*timeStep/60+wellTime
#   })
#   well
# })
# Dot plot with an overlayed box plot -------------------------------------
give.n <- function(x){ ##used to add the N below each sample
  return(c(y = -.15, label = length(x)))
}
# p <- ggplot(data = df, aes(x = factor(strain), y = FL2.A, fill=treatment, ymin=min(FL2.A)*.8, ymax = max(FL2.a)*1.05)) #set up data
# p <- p + geom_point(alpha=0.5, position=position_jitterdodge(dodge.width = 0.9)) # plot transparent (alpha<1) points for each datum
# p <- p + geom_boxplot(outlier.size = 0, position = position_dodge(width = .9)) #add a boxplot overtop to summarize the inner quartiles of the data, leave off outliers as they are plotted above
# p <- p + stat_summary(fun.data = give.n, geom = "text", position=position_dodge(width = .9), size=3, angle=90) #add N below each sample, this could also be modified to include other sample info
# p #show me the plot!


# Steady state analysis ---------------------------------------------------
#' Find significant differences in steady state fluorescence
#' @description Compares several replicates of strain-treatment combinations for statistically significant differences in steady-state fluorescence
#' @param yourFlowSets
#' @param ploidy = the gate to subset your flowsets based on the ploidy of you strains and your cytometer.
#' @param only= 'yeast', 'singlets', or 'doublets'
#' @return
#' @export
#' @examples
steadyState <- function(flowset,ploidy="diploid", only="singlets"){

  ### Number of cells (experiments) in the flowSet
  n_experiment <- length(flowset)

  ### Pulling out data for specific channel to be used
  #channel <- flowSet[,'FL2.A',drop=FALSE]

  ### Gate the samples
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
  ### Convert flowSet to dataframe containing all events for each subset => for plotting
  if (only==F) {

    print("Converting all yeast events...")
    yeastdF <- ddply(pData(yeast), colnames(pData(yeast))[-1],
                     function(tube){
                       fsApply(x = yeast[tube$name],rbind,use.exprs = T)})
    print("Converting doublets events...")
    doubletsdF <- ddply(pData(doublets), colnames(pData(doublets))[-1],
                        function(tube){
                          fsApply(x = doublets[tube$name],rbind,use.exprs = T)})

    print("Converting singlets events...")
    singletsdF<- ddply(pData(doublets), colnames(pData(doublets))[-1],
                       function(tube){
                         fsApply(x = doublets[tube$name],rbind,use.exprs = T)})}
  if (only=="singlets") {
    print("Converting singlets events...")
    singletsdF<- ddply(pData(singlets), colnames(pData(singlets))[-1],
                       function(tube){
                         fsApply(x = singlets[tube$name],rbind,use.exprs = T)})
    return(singletsdF)
  } else if (only=="doublets") {
    print("Converting doublets events...")
    doubletsdF <- ddply(pData(doublets),
                        colnames(pData(doublets))[-1],
                        function(tube){
                          fsApply(x = doublets[tube$name],rbind,use.exprs = T)})
    return(doubletsdF)
  } else if (only=="yeast") {
    print("Converting all yeast events...")
    yeastdF <- ddply(pData(yeast),
                     colnames(pData(yeast))[-1],
                     function(tube){
                       fsApply(x = yeast[tube$name],rbind,use.exprs = T)})
    return(yeastdF)
  }
}

# #########################
# ###  Cytometer Gates  ###
# #########################
#
#
#
# ###These stricter gates were set for the NemLab cytometer "Special Snowflake" 7sept2015 by EPJ
#
# ### EPJyeastGate
# ### Defines an SSC.A vs FSC.A gate.  Includes only the yeast population from a flowSet
#
# yeastGate <<- polygonGate(filterId="Yeast",
#                              .gate=matrix(c(
#                                #xvalues
#                                4e5,50000,40000,1250000,2000000,
#                                #yvalues
#                                10000,15000,0.8e5,6e5,5e5),
#                                ncol=2,nrow=5,dimnames=list(c("1","1","1","1","1"),c("FSC.A","SSC.A"))))
#
# ### Diploid Gates
# ### Diploids are slightly larger and have better separation between singlets/doublets
#
# dipsingletGate <<- polygonGate(filterId="DipSingletGate",
#                                   .gate=matrix(c(
#                                     #x values
#                                     8e4,7e4,8e5,1.25e6,1e6,
#                                     #y values
#                                     8e4,1e5,1e6,1.25e6,1e6),
#                                     ncol=2,nrow=5,dimnames=list(rep(NA,5),c("FSC.A","FSC.H"))))
#
# dipdoubletGate <<- polygonGate(filterId="DipDoubletGate",
#                                   .gate=matrix(c(
#                                     #x values
#                                     4e5,4.5e4,2e5,1.3e6,2e6,2e6,
#                                     #y values
#                                     2.5e5,5e4,2e5,1.3e6,1.5e6,1e6),
#                                     ncol=2,nrow=6,dimnames=list(rep(NA,6),c("FSC.A","FSC.H"))))
#
# ### Haploid Gates -->>still  need to set (next time I process haploid data)
#
# hapsingletGate <<- polygonGate(filterId="HaploidSingletGate",
#                                   .gate=matrix(c(
#                                     #x values
#                                     1e4,6.5e4,15e5,12e5,4e4,
#                                     #y values
#                                     4e4,8e4,16e5,18e5,4e5),
#                                     ncol=2,nrow=5,dimnames=list(rep(NA,5),c("FSC.A","FSC.H"))))
#
# hapdoubletGate <<- polygonGate(filterId="HaploidDoubletGate",
#                                   .gate=matrix(c(
#                                     #x values
#                                     1e4,2e5,30e5,25e5,25e5,6.75e5,
#                                     #y values
#                                     2e4,5e4,20e5,30e5,5e5,8.5e5),
#                                     ncol=2,nrow=6,dimnames=list(rep(NA,6),c("FSC.A","FSC.H"))))

# save(yeastGate, dipsingletGate, dipdoubletGate, hapsingletGate, hapdoubletGate, file = 'data/EPJGates.RData')

# # Dynamic analysis --------------------------------------------------------
# experiment<-"20150605_TIR1AFB2"
# setwd( "/Volumes/Nemlings/Clay/Cytometer/20150605_TIR1AFB2/")
#
# plate0<-read.flowSet(path=paste(experiment,"/",sep=""),alter.names=T)
# sampleNames(plate0)<-paste("0_",sampleNames(plate0),sep="")
# plate1<-read.flowSet(path=paste(experiment,"_1/",sep=""),alter.names=T)
# sampleNames(plate1)<-paste("1_",sampleNames(plate1),sep="")
# testSet<-rbind2(plate0,plate1)
#
# testSet[1]
# pData(testSet)
# head(sampleNames(testSet))
# varLabels(testSet)
# file <- fsApply(testSet,function(x) strsplit(keyword(x)$GUID,".fcs")[[1]])
# # go to file 20150605_TIR1AFB2.R and run sample and treatment assembly
# pData(testSet)$sample<-sample
# pData(testSet)$treatment<-treatment
# phenoData(testSet)
# head(exprs(testSet[[1,'FL2.A']]))
#
# yeast <- Subset(testSet,NemyeastGate)
# singlets <- Subset(yeast,NemdipsingletGate)
# pData(singlets)$name<-rownames(pData(singlets))
# reg<-nplr(exprs(singlets[[1,'Time']]), log(exprs(singlets[[1,'FL2.A']])))
# D(reg)
# plot(exprs(singlets[[1,'Time']]), log(exprs(singlets[[1,'FL2.A']])))
# plot(reg)
# #odd bump late in collection
# #nplr looks good, but now need to apply over each sample:treatment pair to get fits for the degradation
# #nplr also does not like values outside of [0,1] so it would be ideal to normalize
# #perhaps just FL2.A/FSC.A would work
# rownames(pData(singlets))
# colnames(pData(singlets))
# exprs(singlets[["1_G01.fcs",'Time']])
# test<-singlets$"1_A04.fcs"
# head(test$Time)
# head(exprs(test)[,'Time'])
# pData(singlets)[pData(singlets)$name==identifier(test),"time"]
# rep(pData(singlets)[pData(singlets)$name==identifier(test),"time"],length(exprs(test$Time)))
# exprs(test)[,'Time']<-exprs(test)[,'Time']/1000+rep(pData(singlets)[pData(singlets)$name==identifier(test),"time"],length(exprs(test$Time)))
# exprs(test) <- cbind(exprs(test),rep(identifier(test),length(exprs(test$Time))))
# identifier(test)
# num<-0
#
# btime_raw <- fsApply(singlets,function(x)as.numeric(unlist(strsplit(keyword(x)$`$BTIM`,split=":"))))
# btime <-  apply(btime_raw,1,function(x)x[1]*60+x[2]+x[3]/60+x[4]/6000)
# time <- btime-min(btime)
# pData(singlets)$time <- time
#
# pData(singlets)
#
# #T0=F7
# pData(singlets)$time<-pData(singlets)$time-(pData(singlets)[pData(singlets)$name=="0_F07.fcs","time"])
#
# dlply(pData(singlets),.(sample,treatment), function(tube){
#   #subset flowset by name
#   ##first find wells with sample
#   #wells<-intersect(sampleNames(singlets),tube$name)
#   fsApply(x = singlets[tube$name,], function(x){
#     # Get time of each frame in minutes of the day
#     exprs(x)[,'Time']<-exprs(x)[,'Time']/1000+rep(pData(singlets)[pData(singlets)$name==identifier(x),"time"],length(exprs(x$Time)))
#     #exprs(x) <- cbind(exprs(x),'ID'=rep(identifier(x),length(exprs(x$Time))))
#
#   })
#   matrix<-fsApply(x = singlets[tube$name,],rbind,use.exprs = T)
#   #matrix<-fsApply(x = Subset(singlets,as.list(rownames(tube))), FUN = function(x) expr(x))
#   reg<-nplr(matrix[,'Time'], matrix[,'FL2.A']/matrix[,'FSC.A'])
#   reg
#   #num<-num+1
#   #cbind(tube$name, rep(num,length(tube)))
# })
#
# #perhaps instead I should combine each frame from the sample,treatment set into a single frame (I will still need to reset the time for each frame) and then use fsApply to do the nplr
#
# #ddply(pData(singlets),.(sample,treatment),
# #fl_mean <- fsApply(flowset,function(x)mean(x[,channel]),use.exprs=T)
#
#
# #' Linear rate of change of fluorescence
# #' @description  Calculates the rate of degradation or accumulation of fluorescence within the linear range for each strain and treatment combination. Uses nplr package to fit an n-parameter logistic regression curve to each dataset, then computes the second derivative of this curve and finds the range where the second derivative is close to zero. Finally the average of the first derivative across this range is calculated
# #'
# #' @param yourFlowSet, self explanatory.
# #' @param annotation_df, A data frame with all of the wells in the flowset labeled with the strain and treatment in that well.
# #' @param ploidy, the gate to subset your flowsets based on the ploidy of you strains and your cytometer.
# #' @param transform, should the data be transformed and how
# #' @param FLX, which fluorescent channel to use.
# #'
# #' @return df, a dataframe with each strain and treatment and the corresponding mean and standard deviation of the linear degradation rate.
# #' @export
# #'
# #' @examples
# linearRate <- function(yourFlowSet, annotation_df, ploidy, transform=F, only=F, FLX='FL2.A') {
#   # Number of experiments
#   n_experiments <- length(flowset)
#
#
# #Transform the dataset
#   # If using channel="FSC.A", don't use fscanorm
#   if (channel=="FSC.A"&transform=="fscanorm") {
#     print("Channel FSC.A selected with no transform= setting set.")
#     print("Defaulting to no transform (set transform=\"log\" for log transform)")
#     transform=F
#   }
#
#   # Transform FL1.A
#   if (transform != F) {
#     print(paste("Transforming FL1.A using", transform, "transform..."))
#     flowset <- fl1transform(flowset,transform=transform)
#   }
# #Subset the dataset
#   if (ploidy=="haploid") {
#     print("Gating with haploid gates...")
#     yeast <- Subset(flowset,NemyeastGate)
#     singlets <- Subset(yeast,NemhapsingletGate)
#     doublets <- Subset(yeast,NemhapdoubletGate)
#   } else if (ploidy=="diploid") {
#     print("Gating with diploid gates...")
#     yeast <- Subset(flowset,NemyeastGate)
#     singlets <- Subset(yeast,NemdipsingletGate)
#     doublets <- Subset(yeast,NemdipdoubletGate)
#   } else {
#     stop('Error: You must define ploidy="haploid" or ploidy="diploid"')
#   }
#
#   if (only==F) {
#     # Normalize and summarize each subset
#     print("Summarizing all yeast events...")
#     yeastsum <- flsummary(yeast,channel=channel,moments=moments,split=split,transform=transform)
#
#     print("Summarizing doublets events...")
#     doubletsum <- flsummary(doublets,channel=channel,moments=moments,split=split,transform=transform)
#
#     print("Summarizing singlets events...")
#     singletsum <- flsummary(singlets,channel=channel,moments=moments,split=split,transform=transform)
#   } else {
#     if (only=="singlets") {
#       print("Summarizing singlets events...")
#       singletsum <- flsummary(singlets,channel=channel,moments=moments,split=split,transform=transform)
#       return(singletsum)
#     } else if (only=="doublets") {
#       print("Summarizing doublets events...")
#       doubletsum <- flsummary(doublets,channel=channel,moments=moments,split=split,transform=transform)
#       return(doubletsum)
#     } else if (only=="yeast") {
#       print("Summarizing all yeast events...")
#       yeastsum <- flsummary(yeast,channel=channel,moments=moments,split=split,transform=transform)
#       return(yeastsum)
#     } else {
#       print("'only' must be 'singlets','doublets', or 'yeast'")
#       stop()
#     }
#   }
#
#   # Get time of each frame in minutes of the day
#   btime_raw <- fsApply(flowset,function(x)as.numeric(unlist(strsplit(keyword(x)$`$BTIM`,split=":"))))
#   btime <- apply(btime_raw,1,function(x)x[1]*60+x[2]+x[3]/60+x[4]/6000)
#   time <- btime-min(btime)
#
#   #calculate
#   fit <- nplr()
# }
