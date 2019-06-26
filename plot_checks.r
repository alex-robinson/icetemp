library(myr)



# Load data 
if (TRUE) {

    k15a = read.table("data/Kleiner2015_FIg2_EXPA-IIIa-melt.txt",header=TRUE)



    dat = my.read.nc("test.nc")
    dat$time = dat$time*1e-3     # [a] => [ka]

}

fldr  = "plots"
ptype = "png"


# Plot comparison with Kleiner et al. (2015), Exp A
if (TRUE) {

    kb = 1 

    xlim = c(0,300)

    x.at = seq(0,300,by=10)

    myfigure(fldr,"k15expa",asp=0.9,pointsize=14,type=ptype)
    par(xaxs="i",yaxs="i",cex=1.1,cex.lab=1.1,cex.axis=1.5)

    par(mfrow=c(3,1))

    # Panel 1: Tb #############################
    par(plt=c(0.12,0.95,0.02,0.95))
    ylim = c(-30,0)

    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    y.at = pretty(ylim,10)
    abline(h=y.at,v=x.at,col="grey90",lty=1,lwd=0.5)
    axis(2,at=y.at)
    mtext(side=2,line=3.0,las=0,"Basal temp. (C)")

    # Kleiner steady-state lines 
    abline(h=c(-10),col=1,lty=2)
    lines(dat$time,dat$T_pmp[kb,],col=1,lty=2)

    
    lines(dat$time,dat$T_ice[kb,],col=1,lwd=3)

    box() 

    # Panel 2: ab #############################
    par(plt=c(0.12,0.95,0.02,0.95))
    ylim = c(-2.6,3.6)

    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    y.at = pretty(ylim,5)
    abline(h=y.at,v=x.at,col="grey90",lty=1,lwd=0.5)
    axis(2,at=y.at)
    mtext(side=2,line=3.0,las=0,"Melt rate (mm a-1 w.e.)")

    lines(dat$time,-dat$bmb*1e3,col=1,lwd=3)
    
    # Kleiner et al. (2015) analytical solution from Fig. 2
    lines(k15a$time,k15a$ab,col=2,lwd=3)

    box() 
    
    # Panel 2: H_w #############################
    par(plt=c(0.12,0.95,0.20,0.95))
    ylim = c(0,160)

    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    y.at = pretty(ylim,10)
    abline(h=y.at,v=x.at,col="grey90",lty=1,lwd=0.5)
    axis(2,at=y.at)
    mtext(side=2,line=3.0,las=0,"Water thickness (m)")

    # Kleiner peak water value 
    abline(h=c(130),col=1,lty=2)
    
    lines(dat$time,dat$H_w,col=1,lwd=3)
        
    axis(1,at=seq(0,300,by=50))
    mtext(side=1,line=2.1,las=0,"Time (ka)")
    box() 
    


    graphics.off()

}