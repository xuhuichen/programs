plot_421_sed_2001 <- function(l2f){
 # find the observational data
 sedDir <- set_sed_data_dir()
 par(cex=0.8)
 # low states
 tmp = read.table(file=paste(sedDir,"/x_newd1.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 erru = lg_nln+tmp[3]
 errd = lg_nln-tmp[4]
 par(fg="blue")
 errbar(t(lg_nu),t(lg_nln),t(erru),t(errd),add=TRUE,pch=16,col="blue")
 tmp = read.table(file=paste(sedDir,"/g_newd1.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 erru = log10(10**tmp[2]+tmp[3])+l2f
 errd = log10(10**tmp[2]-tmp[3])+l2f
 errbar(t(lg_nu),t(lg_nln),t(erru),t(errd),add=TRUE,pch=16,col="blue")
 # high states
 par(fg="cyan")
 # X-ray
 tmp = read.table(file=paste(sedDir,"/x_newa1.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 erru = lg_nln+tmp[3]
 errd = lg_nln-tmp[4]
 errbar(t(lg_nu),t(lg_nln),t(erru),t(errd),add=TRUE,pch=16,col="cyan")
 tmp = read.table(file=paste(sedDir,"/g_newa1.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 erru = log10(10**tmp[2]+tmp[3])+l2f
 errd = log10(10**tmp[2]-tmp[3])+l2f
 errbar(t(lg_nu),t(lg_nln),t(erru),t(errd),add=TRUE,pch=16,col="cyan")
 par(fg="black")
 par(cex=1.2)
}

 #==========================================
plot_421_sed_mw <- function(l2f){
 # find the observational data
 sedDir <- set_sed_data_dir()
 # Optical 48" for March 2001
 lg_nu = c(14.73, 14.73)
 lg_nln = c(-10.02, -9.82)
 err = 1.45e-11
 errd = log10(10**lg_nln - err)+l2f
 erru = log10(10**lg_nln + err)+l2f
 lg_nln = lg_nln+l2f
 par(fg="grey33")
 errbar(t(lg_nu),t(lg_nln),t(erru),t(errd),add=TRUE,pch=5,col="grey33")
 par(fg="black")
 # BeppoSAX 1998 and 2000 data
 tmp = read.table(file=paste(sedDir,"/sax_98_and_00.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 par(cex=0.7)
 points(t(lg_nu), t(lg_nln), pch=16, col="grey66")
 lg_nln = tmp[3]+l2f
 points(t(lg_nu), t(lg_nln), pch=16, col="grey66")
 lg_nln = tmp[4]+l2f
 points(t(lg_nu), t(lg_nln), pch=16, col="grey66")
 # BeppoSAX 1998 and 2000 data
 tmp = read.table(file=paste(sedDir,"/rxte_01_low_and_high.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 points(t(lg_nu), t(lg_nln), pch=16, col="grey33")
 lg_nln = tmp[3]+l2f
 points(t(lg_nu), t(lg_nln), pch=16, col="grey33")
 lg_nln = tmp[4]+l2f
 points(t(lg_nu), t(lg_nln), pch=16, col="grey33")
 par(cex=1.0)
 # Wise 2010 data
 tmp = read.table(file=paste(sedDir,"/mrk421_wise_high.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 points(t(lg_nu), t(lg_nln), pch=5, col="red")
 # Fermi 18 months data
 tmp = read.table(file=paste(sedDir,"/mrk421_fermi_18m.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1] + 23.384 # From GeV to Hz
 lg_nln = tmp[2]+l2f
 points(t(lg_nu), t(lg_nln), pch=5, col="orange")
 par(cex=1.2)
}

 #====================================================

plot_279_sed <- function(l2f){
   # tmpFun <- get(paste("read_sed_",whichSim,sep=""))
   # tmpFun()
 # find the observational data
 sedDir <- set_sed_data_dir()
 # low state data
 tmp = read.table(file=paste(sedDir,"/3c279_low.dat",sep=""),colClasses=c("numeric"))
 lg_nu = log10(tmp[1])
 lg_nln = log10(tmp[2])+l2f
 # plot the data points
 points(t(lg_nu), t(lg_nln), pch=16, col="blue")

 # high state data, the data is in Luminosity units already
 tmp = read.table(file=paste(sedDir,"/3c279_flare1.dat",sep=""),colClasses=c("numeric"))
 lg_nu = log10(tmp[1])
 lg_nln = log10(tmp[2])
 # plot the data points
 points(t(lg_nu), t(lg_nln), pch=16, col="red")
}

 #==================================================

plot_421_sed_2009 <- function(l2f){
 # find the observational data
 sedDir <- set_sed_data_dir()
 # radio state data
 tmp = read.table(file=paste(sedDir,"/mrk421_radio_2009.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 points(t(lg_nu), t(lg_nln), pch=16, col="gray")
 # UV and optical data
 tmp = read.table(file=paste(sedDir,"/mrk421_uvot_2009.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 points(t(lg_nu), t(lg_nln), pch=16, col="gray")
 # Wise data, not in 2009, but in 2010
 tmp = read.table(file=paste(sedDir,"/mrk421_wise_high.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 points(t(lg_nu), t(lg_nln), pch=16, col="red")
 # Fermi data
 tmp = read.table(file=paste(sedDir,"/mrk421_fermi_2009.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 points(t(lg_nu), t(lg_nln), pch=16, col="red")
 # X-ray data
 tmp = read.table(file=paste(sedDir,"/mrk421_xrt_2009.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 points(t(lg_nu), t(lg_nln), pch=16, col="blue")
 tmp = read.table(file=paste(sedDir,"/mrk421_bat_2009.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 points(t(lg_nu), t(lg_nln), pch=16, col="blue")
 # Magic data
 tmp = read.table(file=paste(sedDir,"/mrk421_magic_2009.dat",sep=""),colClasses=c("numeric"))
 lg_nu = tmp[1]
 lg_nln = tmp[2]+l2f
 points(t(lg_nu), t(lg_nln), pch=16, col="blue")
}

plot_421_lc_ic <- function(dataTimeShift){
 # find the observational data
 lcDir <- set_lc_data_dir()
 # TeV light curve 
 tmp = read.table(file=paste(lcDir,"/SAVE_x_and_tev_tev01_wrbr_hfreb2_pca_5_23_32.dat",sep=""),colClasses=c("numeric"))
 xtt = tmp[1]
 rate = tmp[5]
 err = tmp[6]
 time = xtt/1000 # converting from s to ks
 time = time + dataTimeShift
 rateMax = max(rate)
 rateNorm = rate/rateMax
 points(t(time), t(rateNorm),type="b", lty=2, pch=1, col="blue4",cex=0.5)
}

plot_421_lc_sync <- function(dataTimeShift){
 # find the observational data
 lcDir <- set_lc_data_dir()
 tmp = scan(file=paste(lcDir,"/mrk421_24kev.lc",sep=""),what=character(),skip=37,nmax=2)
 t0d_int = as.double(tmp[2])
 tmp = scan(file=paste(lcDir,"/mrk421_24kev.lc",sep=""),what=character(),skip=38,nmax=2)
 t0d_dec = as.double(tmp[2])
 t0d = t0d_int + trunc(t0d_dec/0.8)
 t0s = ((t0d_int-t0d)+t0d_dec)*86400.0
 # X-ray light curve 2-4 keV
 tmp = read.table(file=paste(lcDir,"/mrk421_24kev.lc",sep=""),colClasses=c("numeric"))
 inTime = tmp[2]
 dt = tmp[3]
 rate = tmp[4]
 err = tmp[5]
 inTime = inTime + t0s
 time = inTime/1000 # converting from s to ks
 time = time + dataTimeShift
 rateMax = max(rate)
 rateNorm = rate/rateMax
 points(t(time), t(rateNorm),type="b", lty=2, pch=1, col="gray40",cex=0.5)
}

plot_279_lc_ic <- function(dataTimeShift){
 # find the observational data
 lcDir <- set_lc_data_dir()
 # Fermi light curve
 tmp = read.table(file=paste(lcDir,"/3c279_fermi.lc",sep=""),colClasses=c("numeric"))
 time = tmp[1]
 lumi = tmp[2]
 time = time + dataTimeShift
 lumiMax = max(lumi)
 lumiNorm = lumi/lumiMax
 points(t(time), t(lumiNorm),type="b", lty=2, pch=1, col="blue4",cex=0.5)
}

plot_279_lc_sync <- function(dataTimeShift){
 # find the observational data
 lcDir <- set_lc_data_dir()
 # Optical light curve H
 tmp = read.table(file=paste(lcDir,"/3c279_H.lc",sep=""),colClasses=c("numeric"))
 time = tmp[1]
 lumi = tmp[2]
 time = time + dataTimeShift
 lumiMax = max(lumi)
 lumiNorm = lumi/lumiMax
 points(t(time), t(lumiNorm),type="b", lty=2, pch=1, col="red",cex=0.5)
 # Optical light curve V
 tmp = read.table(file=paste(lcDir,"/3c279_V.lc",sep=""),colClasses=c("numeric"))
 time = tmp[1]
 lumi = tmp[2]
 time = time + dataTimeShift
 lumiMax = max(lumi)
 lumiNorm = lumi/lumiMax
 points(t(time), t(lumiNorm),type="b", lty=2, pch=1, col="blue",cex=0.5)
}

 #=================================================
set_lc_data_dir <- function(){
 return(paste(Sys.getenv("HOME"),"/lightcurves",sep=""))
}
set_sed_data_dir <- function(){
 return(paste(Sys.getenv("HOME"),"/observ",sep=""))
}
 
