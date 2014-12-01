plot_spt <- function(){
 # setup simulation case
 readSim <- set_read_sim()
 whichSim = readSim$whichSim
 iSed = readSim$iSed
 agnType = readSim$agnType
 l2f = readSim$l2f
 # read data header
 tmp = scan(file=paste(whichSim,"_seds.dat",sep=""),what=character(),skip=0,nmax=5)
 t0 = as.double(tmp[2])
 tmax = as.double(tmp[3])
 dt = as.double(tmp[5])
 tmp = scan(file=paste(whichSim,"_seds.dat",sep=""),what=character(),skip=2,nmax=2)
 sedScaling = as.double(tmp[2])
 cat("Parameters read from data file: \n")
 cat("t0=",t0,"tmax=",tmax,"dt=",dt,"sedScaling=",sedScaling,"\n")
 cat("iSed=",iSed,"\n")

 # read data
 
 tmp <- read.table(file=paste(whichSim,"_seds.dat",sep=""), colClasses=c("numeric"))
 lx_kev = log10(tmp[1])
 lx_ev = lx_kev + 3 # photon energy grid in keV
 lx_hz = lx_kev + 17.384 # photon energy grid in Hz
 flux = tmp[c(2:(length(tmp)-1))]
 flux = flux/sedScaling
 # convert the time for labels
 if(agnType == "FSRQ"){
	secToTime = 86400. # Convert to days
	timeUnit = "days"
 } else {
	secToTime = 1000. # Convert to ks
	timeUnit = "ks"
 }
 timeStart = as.integer(((iSed-1)*dt+t0)/secToTime + 0.5)
 timeEnd = as.integer((iSed*dt+t0)/secToTime + 0.5)
 # Additional Emission Components
 if(agnType == "FSRQ"){
	tmp <- read.table(file="../blackbody.in", colClasses=c("numeric"))
	lxfg_kev = log10(tmp[1])
	fg0 = tmp[2]
	fg <- t(approx(t(lxfg_kev),t(fg0),xout=t(lx_kev),rule=2,ties="ordered")$y)
 cat("\n")
	flux = flux + fg
 }
 # Get nu*F_nu, and 
 lseds = log10(flux) +rep(lx_kev,length(flux))
 # plotting limits
 xlimits = c(11.8,27.6)
 ylimits = c(42.6,46.2)
 if (whichSim == "279"){
 xlimits = c(9.,26.)
 ylimits = c(42.2,48.)
 }
 # plot colors
 color = c("cyan","blue","red","green")
 # set the output eps file
 setEPS()
 postscript(file=paste(whichSim,"_seds.eps",sep=""))
 par(mar=c(4,4,1,3.5)+0.3,cex=1.2)
 # plotting
 plot(t(lx_hz),t(lseds[iSed[2]]),"s",xlim=xlimits,ylim=ylimits, lab=c(6,3,15),col=color[2],xlab=expression(nu(Hz)),ylab=expression(paste(nu,L[nu],"(erg ",s^{-1},")")))
 # plot observational data
 if(whichSim == "421"){
 #plot_421_sed_2009(l2f)
 plot_421_sed_mw(l2f)
 plot_421_sed_2001(l2f)
 } else if (whichSim == "279"){
 plot_279_sed(l2f)
 }
 lines(t(lx_hz),t(lseds[iSed[2]]),"s",col=color[2],lwd=1)
 lines(t(lx_hz),t(lseds[iSed[3]]),"s",col=color[3],lwd=1)
 lines(t(lx_hz),t(lseds[iSed[4]]),"s",col=color[4],lwd=1)
 lines(t(lx_hz),t(lseds[iSed[1]]),"s",col=color[1],lwd=1)
 # add an right axis
 axisTick = seq(trunc(ylimits[1]),trunc(ylimits[2]),1)
 axis(4,at=round(axisTick-l2f)+l2f,labels=round(axisTick-l2f))
 mtext(expression(paste(nu,F[nu],"(erg ",s^{-1},")")),4, line=2.5, cex=1.2)
 #labels
 ySpace = (ylimits[2]-ylimits[1])/20.
 text(xlimits[1],ylimits[2]-0,paste("time=",timeStart[1],"-",timeEnd[1],timeUnit),pos=4, col=color[1])
 text(xlimits[1],ylimits[2]-ySpace,paste("time=",timeStart[2],"-",timeEnd[2],timeUnit),pos=4, col=color[2])
 text(xlimits[1],ylimits[2]-ySpace*2,paste("time=",timeStart[3],"-",timeEnd[3],timeUnit),pos=4, col=color[3])
 text(xlimits[1],ylimits[2]-ySpace*3,paste("time=",timeStart[4],"-",timeEnd[4],timeUnit),pos=4, col=color[4])

 dummy <- dev.off()
}
	

