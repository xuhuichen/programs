movie_es <- function(){

# wait for user input
 default_tmp = 2
 pnth <- readline(paste("enter the power to be used in the electron spectrum: [",default_tmp,"] "))
 pnth <- ifelse(pnth == "",default_tmp,as.integer(pnth))
 default_tmp = 1
 time1 <- readline(paste("First EED file to plot: [",default_tmp,"] "))
 time1 <- ifelse(time1 == "",default_tmp,as.integer(time1))
 default_tmp = 10
 time2 <- readline(paste("Last EED file to plot: [",default_tmp,"] "))
 time2 <- ifelse(time2 == "",default_tmp,as.integer(time2))

# read in the gamma grid file, and the the zone numbers
 x <- read.table("E_e.dat")
 gamma = x + 1
 lg = log10(gamma)
# zone numbers
 tmp = scan(file="E_e.dat",what=character(),skip=1,nmax=5)
 nz = as.integer(tmp[3])
 nr = as.integer(tmp[5])
# zones to skip
 zones = nz*3

# determine the y axis limits
 timestepC = formatC(time2, width=4, flag="0")
 tmp <- read.table(file=paste("smp_",timestepC,".dat",sep=""), colClasses=c(rep("NULL",zones),"numeric") )
 Ng = tmp[1]
 lgNg = log10(Ng*gamma**pnth)
 ymax=max(lgNg)+0.2
 ymin=ymax-2

# Set the postscript enviroment for espe files
 setEPS()
# Loop through many electron files, plot
 for (timestep in time1:time2) {
	timestepC = formatC(timestep, width=4, flag="0")
	# open the eps file
	postscript(file=paste("ele_",timestepC,".eps",sep=""))
	tmp <- read.table(file=paste("smp_",timestepC,".dat",sep=""), colClasses=c(rep("NULL",zones), "numeric") )
	Ng = tmp[1]
	lgNg = log10(Ng*gamma**pnth)
	plot(t(lg),t(lgNg),"l",xlim=c(0,6), ylim=c(ymin,ymax), lab=c(6,3,15), xlab=expression(log(gamma)),
	 ylab= bquote(log(gamma^.(pnth)*N[gamma])))
	#mtext("time=",3,line=1)
	#legend("topleft","time=")
	text(0.1,ymax-0.1,paste("time=",timestep),pos=4)
	# Turn off device driver (to flush output to eps)
	dev.off()
 }
 return()
}



### plot_est ###-------------------------------------
plot_est <- function(){

 # wait for user input
 default_tmp = 2
 pnth <- readline(paste("enter the power to be used in the electron spectrum: [",default_tmp,"] "))
 pnth <- ifelse(pnth == "",default_tmp,as.integer(pnth))
 # read in the gamma grid file
 x <- read.table("E_e.dat")
 gamma = x + 1
 lg = log10(gamma)
# zone numbers
 tmp = scan(file="E_e.dat",what=character(),skip=1,nmax=5)
 nz = as.integer(tmp[3])
 nr = as.integer(tmp[5])
# zones to skip
 zones = nz*3
 # times to plot EED
 timestep = c(80,100,200,300) #c(1,100,900,1700)
 # determine the y axis limits
 timestepC = formatC(timestep[4], width=4, flag="0")
 tmp <- read.table(file=paste("smp_",timestepC,".dat",sep=""), colClasses=c(rep("NULL",zones),rep("numeric",1)) )
 Ng = tmp[1]
 lgNg4 = log10(Ng*gamma**pnth)
 ymax=max(lgNg4)+0.2
 ymin=ymax-2
 # read in the date
 for (i in 1:4) {
	timestepC = formatC(timestep[i], width=4, flag="0")
	tmp <- read.table(file=paste("smp_",timestepC,".dat",sep=""), colClasses=c(rep("NULL",zones),rep("numeric",1)) )
	Ng = tmp[1]
	assign(paste("lgNg",i,sep=""), log10(Ng*gamma**pnth))
 }
 # set the output eps file
 setEPS()
 postscript(file="est.eps")
 par(mar=c(4.1,4,1,3.5)+0.3)
 # plotting
 plot(t(lg),t(lgNg1),"l",xlim=c(0,6), ylim=c(ymin,ymax), lab=c(6,3,15), col="cyan", lty="dotdash", lwd=2, xlab=expression(log(gamma)), ylab= bquote(log(gamma^.(pnth)*N[gamma])))
 lines(t(lg),t(lgNg2),col="green",lty="dashed",lwd=2)
 lines(t(lg),t(lgNg3),col="red",lty="longdash",lwd=3)
 lines(t(lg),t(lgNg4),col="black",lwd=1)
 # labels
 ySpace = 0.08
 text(0,ymax-0,paste("time=",timestep[1]),pos=4, col="cyan")
 text(0,ymax-ySpace,paste("time=",timestep[2]),pos=4, col="green")
 text(0,ymax-ySpace*2,paste("time=",timestep[3]),pos=4, col="red")
 text(0,ymax-ySpace*3,paste("time=",timestep[4]),pos=4, col="black")
 #plot(t(lg),t(lgNg3),"l",xlim=c(0,6), ylim=c(ymin,ymax), lab=c(6,3,15), col.lab="red", xlab=expression(log(gamma)), ylab= bquote(log(gamma^.(pnth)*N[gamma])))
 #text(0.1,ymax-0.1,paste("time=",timestep[3]),pos=3, col="red")
 #plot(t(lg),t(lgNg4),"l",xlim=c(0,6), ylim=c(ymin,ymax), lab=c(6,3,15), col.lab="blue", xlab=expression(log(gamma)), ylab= bquote(log(gamma^.(pnth)*N[gamma])))
 #text(1,ymax-0.1,paste("time=",timestep[4]),pos=3, col="blue")
 dev.off()
}
