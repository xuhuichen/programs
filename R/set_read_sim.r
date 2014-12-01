set_read_sim <- function(){
 default_tmp= "421"
 whichSim <- readline(paste("source name: [",default_tmp,"] "))
 whichSim <- ifelse(whichSim == "",default_tmp,whichSim)
 iSed = c(1,21,28,30)
 # source specific parameters
 # use H=70km/s/Mpc,  (icosmos.co.uk, matter 0.28, dark energy 0.72)
 if(whichSim=="421"){
	iSed=c(11,15,25,29)
	l2f=54.34 # redshift of Mrk421 is 0.0308 (wiki and NED) D_421=0.135Gpc
 } 
 else if(whichSim=="1510"){
	iSed = c(1,9,15,23) 
	l2f=56.65 # redshift of PKS1510 is 0.36 (NED) D_1510=1.93Gpc (2012.1.3 used 56.45, and later 56.42)
 }
 else if(whichSim=="0208"){
	iSed = c(1,5,8,11) 
	l2f=57.73 # redshift of PKS0208 is 1.003 (2008ApJS.175.97) D_0208=6.71Gpc
 }
 else if(whichSim=="279"){
	iSed = c(1,5,8,11) 
	l2f=57.06 # redshift of 3C279 is 0.536 (2010Nature.463.919) D_279=3.10Gpc
 }
 else if(whichSim=="1424"){
	l2f=57.19 # redshift of PKS1424 is 0.6035 (Cerruti 1424 draft) D_1424=3.58Gpc
 }


 if (whichSim=="1510" || whichSim=="279" || whichSim=="0208"){
	agnType = "FSRQ"
 } else {
	agnType = "HBL"
 }
 readSim <- list("whichSim" = whichSim, "iSed"= iSed,"agnType"=agnType,"l2f"=l2f) 
 return(readSim)
}
