## Plot the bubble plot for drought clusters
## 
## Oldrich Rakovec oldrich.rakovec@ufz.de
## Nov. 2019
## ##

## #########################################################################
rm(list = ls(all = TRUE))
Sys.setenv(TZ="UTC") #
graphics.off()

library(zoo)
library(xts)
library(wordcloud)

maindir="/Users/rakovec/Nextcloud/Cloud/UFZ/projects/xeros/joint_of_2018_2019/results/clustering/"
setwd(maindir)

## datadir="/data/xeros/joint_of_2018_2019/mhm_output/merged_1766_2020"
datadir=maindir

## CASES=c("ClTh3_nCells10","ClTh3_nCells30","ClTh3_nCells50","ClTh3_nCells70",
##         "ClTh5_nCells10","ClTh5_nCells30","ClTh5_nCells50","ClTh5_nCells70",
##         "ClTh10_nCells10","ClTh10_nCells30","ClTh10_nCells50","ClTh10_nCells70",
##         "ClTh15_nCells10","ClTh15_nCells30","ClTh15_nCells50","ClTh15_nCells70")

CASES=c("ClTh10_nCells60")


## ####################################
cexpt=1.5
plotName=paste("bubble_plot_v2.pdf",sep="")
pdf(plotName, width = 5, height = 5, family="Helvetica")
par(oma=c(0.25,0.25,0.25,0.25), mar=c(3.5,3.75,0.1,0.1), las=1,  cex.lab=cexpt)#, lwd=2,cex.axis=2.7)
par(mgp = c(2.3, 0.3, 0))
## define plotting regions
margins=matrix(NA,nrow=6,ncol=4)
margins[1,]=c(0.2,0.99,0.3,0.99)
margins[2,]=c(0.2,0.95,0.0,0.25)


for (cas in CASES){

    ## read the data from:
    tabledir=paste(datadir,"/smi_calc_",cas,sep="")
    ## tabledir="/Users/rakovec/Nextcloud/Cloud/UFZ/projects/xeros/joint_of_2018_2019/results/clustering/smi_calc_ClTh10_nCells60"
    all_full=read.table(paste(tabledir,"/results_ADM.txt", sep=""), header=TRUE, na = "-9999")

    all_full$aDD[all_full$aDD<3]=NA
    all_full$aDA[all_full$aDA<0.05]=NA

    all_full = all_full[is.na(all_full$aDA)==FALSE,]
    all_full = all_full[is.na(all_full$aDD)==FALSE,]

    all = all_full

    aDA=all$aDA*100
    aDD=all$aDD
    TDM=all$TDM

    mStart=all$mStart
    mEnd=all$mEnd

    pp=c(1,2,3.5,4.5)
    
    legT = c(1,1,1,1)
    

    ## make the point size corresponding to the area:
    ## pchsize =  sqrt(TDM/aDD)/4
    ## pchsize =  TDM/aDD/50
    pchsize =  TDM/500
    
    print("++++++++++++++++++++++")
    print(cas)
    print(length(which(TDM>3000)))

    ## ll = unique(c(which(aDA>28),which(aDD>8.5)))

    ll = unique(c(which(TDM/aDD>100)))
    ll=ll[c(63,61,56,46,33,24,13)]

    color_transparent <- rep(adjustcolor("#878787", alpha.f = 0.5),length(pchsize))

    colx=rev(c("#377eb8","#4daf4a","#984ea3","#ff7f00","#e31a1c","#a65628","#1a1a1a"))
    
    color_transparent[ll]=colx
    calYM = as.yearmon(1766 + seq(0, 3060)/12)

    calYMb=format(calYM, "%Y")
    calYMb

    print(paste(calYMb[mStart[ll]], TDM[ll]/aDD[ll] ,sep=" xxx  "))   

    plot(aDA,aDD, cex = pchsize, pch = 19, xaxt="n", yaxt="n",
         ylim=c(3,1.1*max(aDD,na.rm=T)),xlim=c(5,1.1*max(aDA,na.rm=T)),
         col=color_transparent, main="", #cas,
         ylab="Mean duration [month]", xlab="Mean area [% of all domain]")
    textplot(aDA[ll], aDD[ll], paste(calYMb[mStart[ll]],"-",substr(calYMb[mEnd[ll]],3,4), sep=""), new = FALSE, show.lines = TRUE,  col=c("#878787",rep("#1a1a1a",6)), cex=1.25)

    axis(1, labels=TRUE, tck=0.025, cex.axis=cexpt)
    axis(2, labels=TRUE, tck=0.025, cex.axis=cexpt)
    axis(3, labels=FALSE, tck=0.025)
    axis(4, labels=FALSE, tck=0.025)

    ## #################
    ## add legend to the figure

    pch0=c(7,3.5,0.5)
    corp=c(7.55,13.5)
    
    points(corp[1],corp[2], col=1, cex=pch0[1])
    points(corp[1],corp[2], col=1, cex=pch0[2])
    points(corp[1],corp[2], col=1, cex=pch0[3], pch=19)

    ## TDM is
    ## pchsize = TDM/500
    TDM0=pch0*500

    text(12.75,14.45, TDM0[1], cex=1)
    text(12.75,13.95, TDM0[2], cex=1)
    text(12.75,13.45, TDM0[3], cex=1)

    segments(x0 = corp[1], y0 = corp[2], x1 = 10.8,     y1 = corp[2], lwd = 1)      
    segments(x0 = corp[1], y0 = corp[2]+0.4, x1 = 10.8,  y1 = corp[2]+0.4, lwd = 1)      
    segments(x0 = corp[1], y0 = corp[2]+0.8, x1 = 10.8,   y1 = corp[2]+0.8, lwd = 1)      
    
    
    ## do 50-year period grouping:
    ## call by G1-G5 (50-year splits)

    iS=1766
    p1=c(iS:(iS+50))
    p2=c((iS+51):(iS+100))
    p3=c((iS+101):(iS+150))
    p4=c((iS+151):(iS+200))
    p5=c((iS+201):(iS+253))

    ii1=which(calYMb[mStart]%in%p1==TRUE)
    ii2=which(calYMb[mStart]%in%p2==TRUE)
    ii3=which(calYMb[mStart]%in%p3==TRUE)
    ii4=which(calYMb[mStart]%in%p4==TRUE)
    ii5=which(calYMb[mStart]%in%p5==TRUE)

    gname = rep(NA, length(aDA))

    gname[ii1] = "G1"
    gname[ii2] = "G2"
    gname[ii3] = "G3"
    gname[ii4] = "G4"
    gname[ii5] = "G5"

    calYMS=calYMb[mStart]
    calYME=calYMb[mEnd]

    ally = cbind(all,calYMS,calYME,gname)

    write.table(ally,file=paste(tabledir,"/results_ADM_yyyymm.txt", sep=""),col.names=T, row.names=F, na = "-9999", quote = FALSE )
    write.table(ally[ll,],file=paste(tabledir,"/results_ADM_yyyymm_big.txt", sep=""),col.names=T, row.names=F, na = "-9999", quote = FALSE )

}

graphics.off()


