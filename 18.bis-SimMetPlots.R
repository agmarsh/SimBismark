
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
print(args)

#-------------------- User Variables --------------------
#rm(list=ls())

PlotTag <- c(args[1])

PlotTag1 <- c("03")
PlotTag2 <- c("0777")
workFolder <- c("002-Bismark/")

D <- read.table(paste(workFolder,PlotTag1,"-",PlotTag2,"/11-BisSimMethylScoreTable-CpG-",PlotTag2,".txt",sep=''),sep='\t',header=T) 
summary(D)


# PLOT 01 --------------------------------------------------------------------------

s <- ggplot(data=D, aes(x=ObsMet, y=ExpMet) ) +
 	 geom_point(aes(colour=ExpMet), size=3, alpha=1) +
 	 scale_x_continuous("Expected %MET") +
 	 scale_y_continuous("Bismark %MET") +
 	 #geom_text(label=D$POS,size=2, hjust=0,vjust=0) +
 	 #stat_smooth(method="lm", colour="red3", fill="grey20", size=0.5, alpha=0.5, formula = y ~ x) +
 	 labs(title=paste(PlotTag2,": Bismark Scoring")) +
 	 scale_colour_gradient("Expected %MET",low="orange", high="blue")
 	 #facet_wrap(~GRP,scales="free_x")
 	 # geom_text(data=eq,aes(x = 25, y = 300,label=V1), parse = TRUE, inherit.aes=FALSE)
s

ggsave(s, file=paste(workFolder,PlotTag1,"-",PlotTag2,"/","Rplot-Exp-Obs-Bismark.png",sep=''), dpi=300)


# PLOT 02 --------------------------------------------------------------------------
u <- ggplot(data=D, aes(x=ObsMet)) +
	geom_histogram(colour="blue4", fill="green3", binwidth= 2) +
	scale_x_continuous(limits=c(0,100)) +
	#geom_vline(x=mErr, colour="blue4", size=2, alpha=0.4) +
	labs(title="Bismark %MET Distribution", x="%MET Calls", y="Count")
u
ggsave(u, file=paste(workFolder,PlotTag1,"-",PlotTag2,"/","Rplot-pMET-distribution.png",sep=''), dpi=300)


# MODEL REGRESSION --------------------------------------------------------------------------

cor(D$ObsMet,D$ExpMet)

# D$Dev <- abs(D$Obs - D$pMET)/(D$pMET)
D$logObs <- log10(D$ObsMet+100)
D$logExp <- log10(D$ExpMet+100)
D$logDev <- abs(D$logObs-D$logExp)/D$logExp
D$Dev <- 10^(100*(D$logDev))

qDev <- quantile(D$Dev, probs=c(0.20, 0.25, 0.5, 0.75, 0.95), na.rm=TRUE)
qDev
mErr = round(qDev[3],2)


u <- ggplot(data=D, aes(x=Dev)) +
	geom_histogram(colour="blue4", fill="green3", binwidth= 0.60) +
	scale_x_continuous(limits=c(0,40)) +
	geom_vline(x=mErr, colour="blue4", size=2, alpha=0.4) +
	geom_text(label=paste("Mean = ",mErr, "%",sep=""), x=mErr+1,y=110,hjust=0,colour="gray40") +
	labs(title="%Error Distribution", x="CpG[i] %Error", y="Error Count")
#u
ggsave(u, file=paste(workFolder,PlotTag1,"-",PlotTag2,"/","Rplot-SimMethyl-ErrorDist.png",sep=''), dpi=300)


#w <- ggplot(data=D, aes(x=ExpMet,y=ObsMet)) +
#	geom_point(aes(colour=Dev), size=3, alpha=0.4) +
#	geom_line(aes(x=pMET,y=pMET), colour="red3", size=2, alpha=0.6) +
#	labs(title=paste("Obs vs Exp %MET: Mean Error Rate = ",mErr,"%",sep=""), x="%MET in Sample DNA", y="GenPro Measured %MET") + 
#	scale_colour_gradient2(name=c("%Error"),limits=c(min(D$Dev),30),low="blue4", mid="yellow3",high="red4") +
#	geom_text(label="(1:1)", x=67, y=75, colour="red4")
#	#stat_smooth(method="lm", colour="red3", fill="grey40", size=1, alpha=0.5, formula = y ~ x) 
##w
#ggsave(w, file=paste(workFolder,PlotTag1,"-",PlotTag2,"/","Rplot-SimMethyl-OBS-vs-EXP.png",dpi=300)



#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------


