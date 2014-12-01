plot_lc <- function() {
 # setup simulation case
 readSim <- set_read_sim()
 whichSim = readSim$whichSim
 iSed = readSim$iSed
 agnType = readSim$agnType
 l2f = readSim$l2f
 # user inputs
 default_tmp = 1
 dataTimeShift <- readline(paste("data time shift: ? [",default_tmp,"] ks "))
 dataTimeShift <- ifelse(dataTimeShift == "",default_tmp,as.numeric(dataTimeShift))
 #==============================
 # data reading and manipulation
 # read data header
 cat("Parameters read from data file: \n")
 eLabelLow = c()
 eRangeLow= c()
 eUnitLow = c()
 eLabelHigh = c()
 eUnitHigh = c()
 eRangeHigh= c()
 iLc = 1
 tmp = scan(file=paste(whichSim,"_lc.dat",sep=""),what=character(),skip=iLc-1,nmax=5)
 while(tmp[1] == "#energy"){
	etmp = as.double(tmp[3])
	if(etmp < 1.e-2){
	 eLabelLow[iLc] = etmp*1e3 
	 eUnitLow[iLc] = "eV"
	}else if (etmp>=1.e-2 && etmp < 1.e3){
 	 eLabelLow[iLc] = etmp
	 eUnitLow[iLc] = "keV"
	}else if(etmp>=1.0e3 && etmp<1.0e5){
 	 eLabelLow[iLc] = etmp*1e-3 
	 eUnitLow[iLc] = "MeV"
	}else if(etmp>=1.0e5 && etmp<1.0e8){
	 eLabelLow[iLc] = etmp*1e-6
	 eUnitLow[iLc] = "GeV"
	}else if(etmp>=1.0e8){
	 eLabelLow[iLc] = etmp*1e-9
	 eUnitLow[iLc] = "TeV"
	}
	eRangeLow[iLc] = etmp

	etmp = as.double(tmp[4])
	if(etmp < 1.e-2){
	 eLabelHigh[iLc] = etmp*1e3 
	 eUnitHigh[iLc] = "eV"
	}else if (etmp>=1.e-2 && etmp < 1.e3){
 	 eLabelHigh[iLc] = etmp
	 eUnitHigh[iLc] = "keV"
	}else if(etmp>=1.0e3 && etmp<1.0e5){
 	 eLabelHigh[iLc] = etmp*1e-3 
	 eUnitHigh[iLc] = "MeV"
	}else if(etmp>=1.0e5 && etmp<1.0e8){
	 eLabelHigh[iLc] = etmp*1e-6
	 eUnitHigh[iLc] = "GeV"
	}else if(etmp>=1.0e8){
	 eLabelHigh[iLc] = etmp*1e-9
	 eUnitHigh[iLc] = "TeV"
	}
	eRangeHigh[iLc] = etmp

	cat("energy range",iLc,": ",eLabelLow[iLc], eUnitLow[iLc], "~", eLabelHigh[iLc], eUnitHigh[iLc],"\n")
	iLc = iLc+1
	tmp = scan(file=paste(whichSim,"_lc.dat",sep=""),what=character(),skip=iLc-1,nmax=5)
 }
 nLc = iLc-1 # number of light curves in the file.
 tmp = scan(file=paste(whichSim,"_lc.dat",sep=""),what=character(),skip=nLc,nmax=4)
 t0 = as.double(tmp[2])
 tmax = as.double(tmp[3])
 dt = as.double(tmp[4])
 tmp = scan(file=paste(whichSim,"_lc.dat",sep=""),what=character(),skip=nLc+2,nmax=2)
 lcScaling = as.double(tmp[2])
 cat("t0=",t0,"tmax=",tmax,"dt=",dt,"lcScaling=",lcScaling,"\n")
 cat("iSed=",iSed,"\n")

 # read data
 tmp <- read.table(file=paste(whichSim,"_lc.dat",sep=""), colClasses=c("numeric"))
 time = tmp[1]
 flux = tmp[c(2:(nLc+1))]
 flux = flux/lcScaling
 count = tmp[c((nLc+2):(2*nLc+1))]
 count = count/lcScaling
 # time conversion
 if(agnType =="FSRQ"){
	secToTime = 86400. # Convert to days
	timeUnit = "days"
 } else {
	secToTime = 1000. # Convert to ks
	timeUnit = "ks"
 }
 time = time/secToTime
 # Additional Emission Components
 if(agnType == "FSRQ"){
	tmp <- read.table(file="../blackbody.in",colClasses=c("numeric"))
	xfg_kev = tmp[1]
	lxfg_kev = log10(xfg_kev)
	fg0 = tmp[2]
	nfg0 = fg0/xfg_kev/(1.602*1e-9)
	lRangeLow = log10(eRangeLow)
	lRangeHigh = log10(eRangeHigh)
	fgLow <- approx(t(lxfg_kev),t(fg0),xout=t(lRangeLow),rule=2,ties="ordered")$y
	print(fgLow)
	fgHigh <- approx(t(lxfg_kev),t(fg0),xout=t(lRangeHigh),rule=2,ties="ordered")$y
	nfgLow <- approx(t(lxfg_kev),t(nfg0),xout=t(lRangeLow),rule=2,ties="ordered")$y
	nfgHigh <- approx(t(lxfg_kev),t(nfg0),xout=t(lRangeHigh),rule=2,ties="ordered")$y
	fgFlux = rep(0,nLc)
	fgCount = rep(0,nLc)
	for(iLc in 1:nLc){
	 ixfg = 2
	 while(xfg_kev[ixfg,] < eRangeHigh[iLc]){
		# integration for the first partial grid of the energy ranges.
		if(xfg_kev[ixfg-1,] < eRangeLow[iLc] && eRangeLow[iLc] < xfg_kev[ixfg,] ){
		 fgFlux[iLc]= fgFlux[iLc] + (xfg_kev[ixfg,]-eRangeLow[iLc])*sqrt(fg0[ixfg,]*fgLow[iLc])
		 fgCount[iLc]= fgCount[iLc] + (xfg_kev[ixfg,]-eRangeLow[iLc])*sqrt(nfg0[ixfg,]*nfgLow[iLc])
		}
		# integration for the main parts of the energy ranges.
		if(eRangeLow[iLc] < xfg_kev[ixfg-1,]){
		 fgFlux[iLc]= fgFlux[iLc] + (xfg_kev[ixfg,]-xfg_kev[ixfg-1,])*sqrt(fg0[ixfg,]*fg0[ixfg-1,])
		 fgCount[iLc]= fgCount[iLc] + (xfg_kev[ixfg,]-xfg_kev[ixfg-1,])*sqrt(nfg0[ixfg,]*nfg0[ixfg-1,])
		}
		ixfg = ixfg+1
	 }
	 # integration for the last partial grid of the energy ranges.
	 fgFlux[iLc]= fgFlux[iLc] + (eRangeHigh[iLc]-xfg_kev[ixfg-1,])*sqrt(fgHigh[iLc]*fg0[ixfg-1,])
	 fgCount[iLc]= fgCount[iLc] + (eRangeHigh[iLc]-xfg_kev[ixfg-1,])*sqrt(nfgHigh[iLc]*nfg0[ixfg-1,])
	flux[iLc] = flux[iLc]+fgFlux[iLc]
	count[iLc] = count[iLc]+fgCount[iLc]
	} # end of for loop iLc
 } # end of if FSRQ
 lFlux = log10(flux)
 lCount = log10(count)
 normFlux = flux
 normCount = count
 fMax = c()
 nMax = c()
 for(iLc in 1:nLc){
	fMax[iLc]= max(flux[iLc])
	nMax[iLc]= max(count[iLc])
	normFlux[iLc]= flux[iLc]/fMax[iLc]
	normCount[iLc] = count[iLc]/nMax[iLc]
 }

 # ==================== 
 # Prepare the plot

 # plot limits
 if(agnType == "FSRQ"){
	timeLimits = c(0,100)
	timeUnits = "days"
 } else {
	timeLimits = c(0,300)
	timeUnits = "ks"
 }
 yLimits = c(0.1,1.5)
 # plot colors
 colors = c("green","red","blue","yellow","cyan","orange","black")
 #=====================
 # begin plotting
 setEPS()
 postscript(file=paste(whichSim,"_lcs.eps",sep=""))
 # Top panel
 par(mar=c(4,4,1,3.5)+0.3,cex=1.2,mfrow=c(2,1))
 par(mar=c(0,4.3,4.3,3.5))
 plot(t(time),t(normFlux[1]),"s", col=colors[1], log="y", xaxt='n', xlim=timeLimits, ylim=yLimits,lab=c(6,3,15),ylab=expression(paste("F/",F[max],"(erg/s/",cm^2,")")) )
 axis(1,col.axis="white")
 # plot observational data
 if(dataTimeShift != 0){
	if(whichSim == "421"){
	 plot_421_lc_ic(dataTimeShift)
	} else if (whichSim == "279"){
	 plot_279_lc_ic(dataTimeShift)
	}
 }
 # plot the 3 light curves in the top panel
 ySpace = (yLimits[2]/yLimits[1])**0.08
 iLc=1
 lines(t(time),t(normFlux[iLc]),"s",col=colors[iLc],lwd=1)
 text(timeLimits[1],yLimits[2]/ySpace,paste(eLabelLow[iLc],eUnitLow[iLc],"-",eLabelHigh[iLc],eUnitHigh[iLc]),pos=4, col=colors[iLc])
 iLc=6
 lines(t(time),t(normCount[iLc]),"s",col=colors[iLc],lwd=1)
 text(timeLimits[1],yLimits[2]/ySpace**2,paste(eLabelLow[iLc],eUnitLow[iLc],"-",eLabelHigh[iLc],eUnitHigh[iLc]),pos=4, col=colors[iLc])
 iLc=7 # For this case only plot when there is possibility of detection.
 if((fMax[iLc]/10^l2f)*(eRangeHigh[iLc]-eRangeLow[iLc])/sqrt(eRangeHigh[iLc]*eRangeLow[iLc])>1.e-17){
	lines(t(time),t(normCount[iLc]),"s",col=colors[iLc],lwd=1)
	text(timeLimits[1],yLimits[2]/ySpace**3,paste(eLabelLow[iLc],eUnitLow[iLc],"-",eLabelHigh[iLc],eUnitHigh[iLc]),pos=4, col=colors[iLc])
 }

 # Bottom panel
 par(mar=c(4.3,4.3,0,3.5))
 plot(t(time),t(normFlux[2]),"s", col=colors[2], log="y", xlim=timeLimits, ylim=yLimits,lab=c(6,3,15),xlab=paste("time(",timeUnits,")",sep=""),ylab=expression(paste("F/",F[max],"(erg/s/",cm^2,")")) )
 iLc=2
 lines(t(time),t(normFlux[iLc]),"s",col=colors[iLc],lwd=1)
 text(timeLimits[1],yLimits[2]/ySpace,paste(eLabelLow[iLc],eUnitLow[iLc],"-",eLabelHigh[iLc],eUnitHigh[iLc]),pos=4, col=colors[iLc])
 iLc=3
 lines(t(time),t(normFlux[iLc]),"s",col=colors[iLc],lwd=1)
 text(timeLimits[1],yLimits[2]/ySpace**2,paste(eLabelLow[iLc],eUnitLow[iLc],"-",eLabelHigh[iLc],eUnitHigh[iLc]),pos=4, col=colors[iLc])

 # plot observational data
 if(dataTimeShift != 0){
	if(whichSim == "421"){
	 plot_421_lc_sync(dataTimeShift)
	} else if (whichSim == "279"){
	 plot_279_lc_sync(dataTimeShift)
	}
 }

 dummy <- dev.off()
}
