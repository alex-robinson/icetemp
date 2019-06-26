library(myr)

rho_ice = 910.0 
rho_w   = 1000.0 
cp      = 2009.0        # [J kg-1 K-1]

# Load data 
if (TRUE) {

    experiment = "k15expb"
    is_celcius = FALSE 

    
    filename = paste0("test_",experiment,".nc")

    dat = my.read.nc(filename)
    dat$time = dat$time*1e-3     # [a] => [ka]

    # Calculate melt rate [m/a w.e.]
    dat$bmb_we = dat$bmb*rho_ice/rho_w 

    if (! is_celcius) {
        # Convert to Celcius for plots 
        T0 = 273.15 
        dat$T_ice   = dat$T_ice   - T0 
        dat$T_pmp   = dat$T_pmp   - T0 
        dat$T_srf   = dat$T_srf   - T0 
        dat$T_robin = dat$T_robin - T0     
    }

    # Read data for comparison 
    k15a = read.table("data/Kleiner2015_FIg2_EXPA-IIIa-melt.txt",header=TRUE)

}

fldr  = "plots"
ptype = "png"


# Plot comparison with Kleiner et al. (2015), Exp A
if (TRUE & experiment == "k15expa") {

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
    mtext(side=2,line=3.0,las=0,"Basal temp. (°C)")

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

    lines(dat$time,-dat$bmb_we*1e3,col=1,lwd=3)
    
    # Kleiner et al. (2015) analytical solution from Fig. 2
    lines(k15a$time,k15a$ab,col=2,lwd=3)

    box() 
    
    # Panel 2: H_w #############################
    par(plt=c(0.12,0.95,0.20,0.95))
    ylim = c(-5,160)

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

# Plot comparison with Kleiner et al. (2015), Exp B
if (TRUE & experiment == "k15expb") {

    kt = length(dat$time) 

    T_ref = 223.15   # [K] 
    E_ref = (T_ref*rho_ice*cp) / rho_ice *1e-3 
    

    ylim = c(-0.001,1)
    y.at = seq(0,1,by=0.1)

    myfigure(fldr,"k15expb",asp=2.3,pointsize=12,type=ptype)
    par(xaxs="i",yaxs="i",cex.axis=1.0)

    par(mfrow=c(1,3))

    # Panel 1: Tb #############################
    par(plt=c(0.12,0.95,0.12,0.95))
    xlim = c(92,108)
    x.at = c(92,96,100,104,108)
    
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    abline(h=y.at,v=x.at,col="grey90",lty=1,lwd=0.5)
    axis(1,at=x.at)
    axis(2,at=y.at)
    mtext(side=1,line=1.6,las=0,cex=0.7,"Enthalpy (kJ kg-1)")

    lines(dat$enth[,kt]/rho_ice*1e-3 - E_ref,dat$zeta,col=1,lwd=3)

    box() 

    # Panel 2: ab #############################
    par(plt=c(0.12,0.95,0.12,0.95))
    xlim = c(-3.2,0.5)
    x.at = pretty(xlim,5)
    
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    abline(h=y.at,v=x.at,col="grey90",lty=1,lwd=0.5)
    axis(1,at=x.at)
    axis(2,at=y.at)
    mtext(side=1,line=1.6,las=0,cex=0.7,"Temperature (°C)")

    lines(dat$T_ice[,kt],dat$zeta,col=1,lwd=3)

    box() 
    
    # Panel 2: H_w #############################
    par(plt=c(0.12,0.95,0.12,0.95))
    xlim = c(-0.1,3)
    x.at = pretty(xlim,10)
    
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    abline(h=y.at,v=x.at,col="grey90",lty=1,lwd=0.5)
    axis(1,at=x.at)
    axis(2,at=y.at)
    mtext(side=1,line=1.6,las=0,cex=0.7,"Water content (%)")

    # Kleiner peak water value 
    abline(h=c(130),col=1,lty=2)
    
    lines(dat$omega[,kt]*100,dat$zeta,col=1,lwd=3)

    box() 
    


    graphics.off()

}