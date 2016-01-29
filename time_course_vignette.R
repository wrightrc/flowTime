##Nemhauser Lab Cytometer Settings
##load in the packages and codes that R needs to analyze cytometer data
library(flowViz)
library(flowCore)
library(ggplot2)

##Set the title of the experiment,
##which should be the name of the folder containing your flow data
experiment<-"time_course_vignette/20130319_OFC9"
setwd(experiment)

#how does this data look?
expt1<-read.csv(paste0(experiment,'.csv'), row.names=1)
qplot(x = time, y= FL1.Amean, data = expt1, colour = strain, fill = factor(treatment)) + geom_line()
#strain 3 looks really nice, I will subset this.

#read in data
plate0<-read.flowSet(path=paste(experiment,"_0/",sep=""),alter.names=T)
sampleNames(plate0)<-paste("0_",sampleNames(plate0),sep="")
plate1<-read.flowSet(path=paste(experiment,"_1/",sep=""),alter.names=T)
sampleNames(plate1)<-paste("1_",sampleNames(plate1),sep="")
plate2<-read.flowSet(path=paste(experiment,"_2/",sep=""),alter.names=T)
sampleNames(plate2)<-paste("2_",sampleNames(plate2),sep="")
plate3<-read.flowSet(path=paste(experiment,"_3/",sep=""),alter.names=T)
sampleNames(plate3)<-paste("3_",sampleNames(plate3),sep="")
plate4<-read.flowSet(path=paste(experiment,"_4/",sep=""),alter.names=T)
sampleNames(plate4)<-paste("4_",sampleNames(plate4),sep="")

dat<-rbind2(plate0,plate1)
dat<-rbind2(dat,plate2)
dat<-rbind2(dat,plate3)
dat<-rbind2(dat,plate4)

write.flowSet(dat,outdir=experiment)
annotation <- expt1[,c("strain","RD","ARF","AFB","treatment")]
annotation$name <- paste0(row.names(annotation),'.fcs')
adat <- annotateFlowSet(dat, annotation_df = annotation)
pData(adat)
adat <- subset(adat, subset = pData(adat)$strain == 3)
write.flowSet(adat, outdir = 'inst/extdata/tc_example')
annotation <- subset(annotation, annotation$strain == 3)
write.csv(annotation, 'inst/extdata/tc_example.csv')

dat_sum<-summary.cyt(adat, ploidy="diploid", only="singlets", channel="FL1.A")
colnames(dat_sum)[9]<-'X'
join(dat_sum, annotation,'X')
##need to reset gates, but these are bound in our environment
ls.str(envir = 'package:flowTime')
exists('yeastGate')
ls('yeastGate')
where('yeastGate')
unlockBinding('yeastGate', 'package:flowTime')
rm('yeastGate')
flowTime:::yeastGate
##I think I will just move these to the extdata folder, that way they won't get loaded automatically.
#I can create a default set that will be loaded prior to any summary function
#I will also write a function to create new gateSets and then they can be set to default


#OK, so I have gone back and edited summary() to handle the smoother annotation functions I've written, and this works.
#Now I need to work on gateCreator
#save the Gates I currently have in memory
saveGates(fileName = 'FCGates.RData')



################
#YPH
################
dat_sum <- YPHsummary.cyt(dat, ploidy = "diploid", only = 'singlets', channel='FL3.A')


dat_sum$FBP = rep(c("TIR1","AFB5","AFB2","TIR1-AFB2","TIR1-AFB5","AFB2-TIR1","AFB2-AFB5","AFB5-TIR1","AFB5-AFB2"),length(dat_sum$FBP)/9)
?rep
dat_sum$treatment = c(rep(c(rep("mock",9),rep("10uM IAA",9)),3),rep("10uM IAA",9),rep("mock",9),rep("10uM IAA",9),rep("10uM IAA",9),rep("10uM IAA",9),rep("10uM IAA",9),rep("mock",9),rep("10uM IAA",9),rep("mock",9),rep("10uM IAA",9))


dat_sum_blue
setwd(experiment)
write.csv(dat_sum,paste(experiment,".csv",sep=""))
#write.csv(dat_sum_blue,"20140411_full_circuit_blue.csv")

##open the .csv file in excel to label,then save it as .csv again and reopen it in R

experiment<-"20140618_F-box_swaps"
dat_sum<-read.csv(paste(experiment,"_T0.csv",sep=""),row.names=1)

dat_sum
dat_sum[1,10]<-'Substrate'

##set the t=0 point, must know which well/timepoint auxin was added
dat_sum[dat_sum$file=="A01","btime"]
dat_sum$time=dat_sum$btime-dat_sum[dat_sum$file=="1B01","btime"]
dat_sum[dat_sum$file=="0_E12","btime"]-3
dat_sum
write.csv(dat_sum,paste(experiment,"_T0.csv", sep=""))

##background subtraction, requires that there be strains labeled "noYFP"
dat_bs<-addbs(dat_sum, column="FL3.Amean")

write.csv(dat_bs,"20140421_COI1_JAZ_t0_bs.csv")
write.csv(dat_bs_blue,"20140411_fullcircuit_bs.csv")

dat_bs_blue<-read.csv("20140411_fullcircuit_bs.csv")

#Bar plot of a single time point
dat_sum<-read.csv(paste(experiment,"_T0.csv",sep=""))
FL3mean<-dat_sum$FL3.Amean
par(las=2) # make label text perpendicular to axis
par(mar=c(8,8,4,2)) # increase y-axis margin.
barplot(height=FL3mean, main="Average Fluorescence", axes=T, horiz=T, names.arg=dat_sum$Strain)
FL3median<-dat_sum$FL3.Amedian
barplot(height=FL3median, main="Median Fluorescence", axes=T, horiz=T, names.arg=dat_sum$Strain)


qplot(data=dat_sum,x=time,y=FL3.Amean_bs, colour=factor(strain),shape=factor(treatment)) + geom_line() + theme_clean()+ geom_point(size=I(5))+facet_grid(~afb)->plot1
plot1


qplot(data=subset(dat_bs_blue, afb!="AFB2" & strain!="mCerulean IAA1.Fl"),x=time,y=FL4.Amean_bs, colour=factor(strain),shape=factor(treatment)) + geom_line() + theme_clean()+ geom_point(size=I(5))+facet_grid(~afb)->plot2
plot2

qplot(data=subset(dat_bs_blue, afb!=""),x=time,y=FL4.Amean_bs, colour=factor(strain),shape=factor(treatment)) + geom_line() + theme_clean()+ geom_point(size=I(5))+facet_grid(~afb)->plot2
plot2


qplot(data=subset(dat_sum),x=time,y=FL3.Amean_bs, colour=factor(strain),shape=factor(treatment)) +geom_line () + theme_clean()+ geom_point(size=I(5))+facet_grid(~afb)->plot1
plot1

qplot(data=subset(dat_bs, afb!="PpAFB4A" & afb!="AFB2" & strain!="PpIAA1a" & strain!="PpIAA1b" & strain!="PpIAA2"),x=time,y=FL3.Amean_bs, colour=factor(strain),shape=factor(treatment)) +geom_line () + theme_clean()+ geom_point(size=I(5))+facet_grid(~afb)->plot1
plot1


qplot(data=subset(dat_bs, afb!="TIR1 C140A" & afb!="TIR1 G142A" & afb!="PpAFB4A" & strain!="IAA7" & strain!="YFP" & strain!="IAA28" & strain!="IAA7 "),x=time,y=FL3.Amean_bs, colour=factor(strain),shape=factor(treatment)) +geom_line () + theme_clean()+ geom_point(size=I(5))+facet_grid(~afb)->plot1
plot1


qplot(data=subset(dat_bs_blue, strain!="iLOV" & strain!= "iLOVKR" & strain!="Cerulean" & strain!="mCer3" & strain!="mTur2"),x=time,y=FL3.Amean_bs, colour=factor(strain),shape=factor(treatment)) + geom_line() + theme_clean()+ geom_point(size=I(5))->plot3
plot3
qplot(data=subset(dat_bs_blue, strain!="iLOV" & strain!= "iLOVKR" & strain!="Cerulean" & strain!="mCer3" & strain!="mTur2"),x=time,y=FL4.Amean_bs, colour=factor(strain),shape=factor(treatment)) + geom_line() + theme_clean()+ geom_point(size=I(5))->plot4
plot4



qplot(data=subset(dat_bs_blue, strain!="YFPcontrol" & strain!= "fullcircuit" & strain!= "Ceruleancontrol"),x=time,y=FL4.Amean_bs, colour=factor(strain),shape=factor(treatment)) + geom_line() + theme_clean()+ geom_point(size=I(5))+facet_grid(~afb)->plot2
plot2


##normalization
##need to know the names of all your label columns
##typically use method=3
dat_norm<-addnorm(dat_bs,factor_in=c("afb","strain","treatment"),method=3,column="FL3.Amean_bs")
dat_norm
write.csv(dat_norm,"20140321_PpIAAsAFBs.csv")


qplot(data=subset(dat_norm, afb!="PpAFB1" & afb!="PpAFB2" & afb!="PpAFB3" & afb!="PpAFB4" & strain!="IAA7" &strain!="IAA28" & strain!="YFP"),x=time,y=FL3.Amean_bs, colour=factor(strain),shape=factor(treatment)) + geom_line() + theme_clean()+ geom_point(size=I(6)) +facet_grid(~afb)->plot2
plot2

qplot(data=subset(dat_bs_blue, laser!="on" ),x=time,y=FL3.Amean_bs, colour=factor(strain),shape=factor(treatment)) + geom_line() + theme_clean()+ geom_point(size=I(6)) +facet_grid(~afb)->plot2
plot2

qplot(data=subset(dat_bs_blue, laser!="on " & strain!="Venus IAA1.Fl" & strain!="YFP IAA1.Fl" & strain!="mTur2 IAA1.Fl" & strain!="Venus"),x=time,y=FL3.Amean_bs, colour=factor(strain),shape=factor(treatment)) + geom_line() + theme_clean()+ geom_point(size=I(6)) +facet_grid(~afb)->plot2
plot2

qplot(data=subset(dat_norm),x=time,y=FL3.Amean_bs, colour=factor(strain),shape=factor(dilution.factor)) + geom_line() + theme_clean()+ geom_point(size=I(6)) +facet_grid(~shaker)->plot1
plot1


qplot(data=subset(dat_bs_blue, laser!="off"),x=time,y=FL4.Amean_bs, colour=factor(strain),shape=factor(treatment)) + geom_line() + theme_clean()+ geom_point(size=I(5))->plot4
plot4


qplot(data=subset(dat_bs_blue),x=time,y=FL3.Amean_bs, colour=factor(strain),shape=factor(treatment)) + geom_line() + theme_clean()+ geom_point(size=I(5))->plot2
plot2





qplot(data=subset(dat_norm,strain!="noYFP" & treatment!="0"),x=time,y=normed,colour=factor(strain),shape=factor(afb),linetype=factor(treatment)) + geom_line() + theme_bw() + facet_wrap(~label) + geom_point(size=I(5))->plot3
plot3
ggsave(plot3,file="",height=10, width=10)

qplot(data=subset(dat_norm,strain!="noYFP" & strain!="YFP" & treatment!="0"),x=time,y=normed,colour=factor(label),shape=factor(strain)) + geom_line() + theme_bw()+ facet_wrap(~strain)+geom_point(size=I(5))->plot4
plot4
ggsave(plot4,file="20131223_plot4.png",height=10, width=10)


qplot(data=subset(dat_norm,strain!="noYFP" & strain!="IAA14.FL" & strain!="IAA14.t1" & strain!="IAA14.t2"),x=time,y=normed,colour=factor(treatment),shape=factor(truncation)) + geom_line() + theme_bw()+ geom_point(size=I(5))->plot4
plot4
plot4
ggsave(plot4,file="20130718_aux1/fold_init.png", height=10, width=20)

qplot(data=subset(dat_norm,strain!="noYFP"),x=time,y=normed, colour=factor(label),shape=factor(treatment)) + geom_line() + theme_bw()+ geom_point(size=I(5))->plot5
plot5
ggsave(plot5,file="20130926_IAA1_trunc_full_rep1/plot5.png", height=10, width=10)

plot2
qplot(data=dat_norm,x=time,y=normed, colour=factor(treatment),shape=factor(afb)) + geom_line() + facet_wrap(~label) + xlim(c(0,120)) + theme_bw()+ geom_point(size=I(5))
