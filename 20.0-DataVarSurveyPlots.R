
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
print(args)

#-------------------- User Variables --------------------
rm(list=ls())

PlotTag <- c(args[1])

PlotTag1 <- c("01")
PlotTag <- c("0602")
workFolder <- c("002-CrunchSeq/")

D <- read.table(paste(workFolder,PlotTag1,"-",PlotTag,"/","06-SAMscore-Rtable.txt",sep=''),sep='\t',header=T) 
summary(D)
D$LR <- ((D$Left - D$Right)**2)**0.5
D$LR <- D$Left/D$Right
D$LR[which(D$LR < 1.0)] <- D$Right[which(D$LR < 1.0)]/D$Left[which(D$LR < 1.0)]

# PLOT 02 --------------------------------------------------------------------------
s <- ggplot(data=D, aes(x=UMTcuts, y=PeakCOV) ) +
	 geom_point(aes(colour=RefMet), size=3, alpha=1) +
	 labs(title=PlotTag, x='RefMet',y='UMTcuts') +
	 scale_colour_gradient("%MET", low="green3",high="orange3")
s

ggsave(s, file="Rplot-RefMet-UMTcuts.png",dpi=300)


# PLOT 01 --------------------------------------------------------------------------

p <- ggplot(data=D,aes(x=RefMet,y=ObsMet)) + 
	 geom_abline(slope=1, intercept=0, colour="gray50", alpha = 0.5) +
	 #geom_line(aes(x=c(0,100),y=c(0,100)), colour="gray50", alpha = 0.5) +
	 geom_point(aes(colour=RefMet), size=3.5, alpha=1) +
	 scale_y_continuous(limits=c(0,100)) +
	 labs(title="MULTIPLE CpG CUT", x='Expected %MET',y='GP %Methylation Score') +
	 scale_colour_gradient2("%MET")
p

ggsave(p, file=paste(workFolder,PlotTag1,"-",PlotTag2,"/","Rplot-ExpObs-pMET-MANYCUT.png",sep=''), dpi=300)



# # PLOT 03--------------------------------------------------------------------------
# t <- ggplot(data=D, aes(x=UMTscore, y=pMET, group=GRP) ) +
 	 # geom_point(aes(colour=GRP), size=5, alpha=0.1) +
 	 # scale_x_continuous("Unmethylated Score Metric") +
 	 # scale_y_continuous("% CpG MET") +
 	 # stat_smooth(method="lm", colour="red3", fill="grey40", size=1, alpha=0.5, formula = y ~ x) +
 	 # labs(title=paste(PlotTag,": UN-methylated Score Metric by GROUP")) +
 	 # facet_wrap(~GRP,scales="free_x")
# #t
# ggsave(t, file="Rplot-SimMethyl-UMTscore-GROUP.png",dpi=300)


# PLOT 04 --------------------------------------------------------------------------
D$adjMET <- D$GPmet
q <- ggplot(data=D, aes(x=adjMET, y=pMET) ) +
	 geom_point(aes(colour=COV),size=3, alpha=0.3) +
	 scale_x_continuous("Methylated Score Metric") +
	 scale_y_continuous("% CpG MET") +
	 labs(title=PlotTag) +
	 scale_colour_gradient(limits=c(min(D$COV),max(D$COV)),low="yellow2",high="red4")
#q
ggsave(q, file="Rplot-SimMethyl-METscore-vs-percentM.png",dpi=300)


# PLOT 05 --------------------------------------------------------------------------
D$adjUMT <- D$GPumt
r <- ggplot(data=D, aes(x=adjUMT, y=pMET) ) +
	 geom_point(aes(colour=COV),size=3, alpha=0.3) +
	scale_x_continuous("Unmethylated Score Metric") +
	scale_y_continuous("% CpG MET") +
	labs(title=PlotTag) +
	scale_colour_gradient(limits=c(min(D$COV),max(D$COV)),low="yellow2",high="red4")
#r
ggsave(r, file="Rplot-SimMethyl-UMTscore-vs-percentM.png",dpi=300)


# PLOT 06 --------------------------------------------------------------------------
v <- ggplot(data=D, aes(x=COV, y=pMET) ) +
 	 geom_point(aes(colour=GRP),size=3, alpha=0.4) +
 	scale_x_continuous("CpG Seq Coverage") +
 	scale_y_continuous("% CpG MET") +
 	labs(title=PlotTag) 
#v
ggsave(v, file="Rplot-SimMethyl-COV-vs-percentM.png",dpi=300)


## PLOT 07 --------------------------------------------------------------------------
#x <- ggplot(data=D, aes(x=adjUMT,y=adjMET)) +
#	geom_point(aes(colour=pMET), size=3, alpha=0.4) +
#	# geom_line(aes(x=pMET,y=pMET), colour="red3", size=2, alpha=0.6) +
#	labs(title="Known Methylation", x="UMTscore", y="METscore") + 
#	scale_colour_gradient2(name=c("Expected %MET"),limits=c(0,100),low="white",high="blue3")
#	#geom_text(label="(1:1)", x=67, y=75, colour="red4")
#	#stat_smooth(method="lm", colour="red3", fill="grey40", size=1, alpha=0.5, formula = y ~ x) 
##x
#ggsave(x, file="Rplot-SimMethyl-adjUMT-vs-adjMET.png",dpi=300)



# MULTIPLE LINEAR REGRESSION --------------------------------------------------------------------------
# MULTIPLE LINEAR REGRESSION --------------------------------------------------------------------------

D$sqPH <- D$PeakHT^2
D$sqPC <- D$PeakCOV^2
D$sqUMT <- D$UMTcuts^2
regmodel <- lm(RefMet ~ PeakHT + PeakCOV + UMTcuts + sqPH + sqPC + sqUMT, data=D)
regmodel
summary(regmodel)
anova(regmodel)
D$Obs <- predict(regmodel)
cor(D$RefMet,D$Obs)
# D$Dev <- abs(D$Obs - D$pMET)/(D$pMET)
D$logObs <- log10(D$Obs+100)
D$logExp <- log10(D$RefMet+100)
D$logDev <- abs(D$logObs-D$logExp)/D$logExp
D$Dev <- 10^(100*(D$logDev))

qDev <- quantile(D$Dev, probs=c(0.20, 0.25, 0.5, 0.75, 0.95), na.rm=TRUE)
qDev
mErr = round(qDev[3],2)
mErr

u <- ggplot(data=D, aes(x=Dev)) +
	geom_histogram(colour="blue4", fill="green3", binwidth= 0.60) +
	scale_x_continuous(limits=c(0,40)) +
	geom_vline(x=mErr, colour="blue4", size=2, alpha=0.4) +
	labs(title="%Error Distribution", x="CpG[i] %Error", y="Error Count")
u <- u +
	geom_text(label=paste("Mean = ",mErr, "%",sep=""), x=mErr+1,y=max(ggplot_build(u)$data[[1]]$count)*0.60,hjust=0,colour="gray40")

u
ggsave(u, file="Rplot-SimMethyl-ErrorDist.png",dpi=300)

D$pMT2 <- D$RefMet^2
regline <- lm(Obs ~ RefMet, data=D )
summary(regline)
D$ObsRegr <- predict(regline)


w <- ggplot(data=D,aes(x=RefMet,y=Obs)) +
	geom_point(aes(colour=Dev), size=3, alpha=0.4) +
	geom_line(aes(x=RefMet,y=RefMet), colour="red3", size=1, alpha=0.6) +
	labs(title=paste("Obs vs Exp %MET: Mean Error Rate = ",mErr,"%",sep=""), x="%MET in Sample DNA", y="GenPro Measured %MET") + 
	scale_colour_gradient2(name=c("%Error"),limits=c(min(D$Dev),20),low="blue4", mid="yellow3",high="red4")
	
w
ggsave(w, file="Rplot-SimMethyl-OBS-vs-EXP.png",dpi=300)


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------



