library(myr)

rho_ice  = 910.0 
rho_w    = 1000.0 
cp       = 2009.0        # [J kg-1 K-1]
L_ice    = 333500.0      # [J kg-1] Latent heat

convert_to_enthalpy = function(T_ice,omega,T_pmp,cp) {
    
    enth = (1.0-omega)*(rho_ice*cp*T_ice) + omega*(rho_w*(cp*T_pmp+L_ice))

    return(enth)
}

# Load data 
if (TRUE) {

    experiment = "k15expb"
    is_celcius = FALSE 

    T_ref = 223.15   # [K] 
    E_ref = (T_ref*cp)*1e-3   # [kJ kg-1]
    
    filename = paste0("test_",experiment,".nc")

    dat = my.read.nc(filename)
    dat$time = dat$time*1e-3     # [a] => [ka]

    # Calculate melt rate [m/a w.e.]
    dat$bmb_we = dat$bmb*rho_ice/rho_w 

    # Adjust enthalpy for reference
    dat$enth = dat$enth*1e-3 - E_ref 

    # Calculate enthalpy offline 
    enth = convert_to_enthalpy(dat$T_ice,dat$omega,dat$T_pmp,dat$cp)
    # enth = enth/rho_ice*1e-3 - E_ref 
    enth = enth*1e-3 - E_ref 
    
    if (! is_celcius) {
        # Convert to Celcius for plots 
        T0 = 273.15 
        dat$T_ice   = dat$T_ice   - T0 
        dat$T_pmp   = dat$T_pmp   - T0 
        dat$T_srf   = dat$T_srf   - T0 
        dat$T_robin = dat$T_robin - T0     
    }

    # Read data for comparison 
    k15a   = read.table("data/Kleiner2015/Kleiner2015_EXPA_Fig2-IIIa-melt.txt",header=TRUE)
    k15b   = read.table("data/Kleiner2015/Kleiner2015_EXPB_analytic_nz401_z.dat",header=TRUE)
    k15b$zeta = k15b$z / max(k15b$z)
    k15bts = read.table("data/Kleiner2015/Kleiner2015_EXPB_extra-series_enth_b.txt",header=TRUE)
    #k15bts$enth = k15bts$enth + (173.15-T_ref)*cp *1e-3  # Correction for T_ref=173.15 in this data 
    k15bts$enth = k15bts$enth + (k15b$enth[1]-max(k15bts$enth))

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

    ylim = c(0,1)
    y.at = seq(0,1,by=0.1)

    col = c("black","#d6604d")
    lwd = c(2,2)
    lty = c(1,2)

    myfigure(fldr,"k15expb",asp=2.3,pointsize=12,type=ptype)
    par(xaxs="i",yaxs="i",cex.axis=1.0)

    par(mfrow=c(1,3))

    # Panel 1: Enthalpy #############################
    par(plt=c(0.12,0.95,0.12,0.95))
    xlim = c(92,108)
    x.at = c(92,96,100,104,108)
    
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    abline(h=y.at,v=x.at,col="grey90",lty=1,lwd=0.5)
    axis(1,at=x.at)
    axis(2,at=y.at)
    mtext(side=1,line=1.6,las=0,cex=0.7,"Enthalpy (kJ kg-1)")

    # Kleiner et al. (2015) analytical solution from Fig. 4
    lines(k15b$enth,k15b$zeta,col=col[1],lwd=lwd[1],lty=lty[1])

    lines(dat$enth[,kt],dat$zeta,col=col[2],lwd=lwd[2],lty=lty[2])

    # Caculated offline...
    lines(enth[,kt],dat$zeta,col=3,lwd=1,lty=1)

    box() 

    # Panel 2: T_ice #############################
    par(plt=c(0.12,0.95,0.12,0.95))
    xlim = c(-3.2,0.5)
    x.at = pretty(xlim,5)
    
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    abline(h=y.at,v=x.at,col="grey90",lty=1,lwd=0.5)
    axis(1,at=x.at)
    axis(2,at=y.at)
    mtext(side=1,line=1.6,las=0,cex=0.7,"Temperature (°C)")

    # Kleiner et al. (2015) analytical solution from Fig. 4
    lines(k15b$T_ice,k15b$zeta,col=col[1],lwd=lwd[1],lty=lty[1])

    lines(dat$T_ice[,kt],dat$zeta,col=col[2],lwd=lwd[2],lty=lty[2])

    box() 
    
    # Panel 2: omega #############################
    par(plt=c(0.12,0.95,0.12,0.95))
    xlim = c(-0.1,3)
    x.at = pretty(xlim,10)
    
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    abline(h=y.at,v=x.at,col="grey90",lty=1,lwd=0.5)
    axis(1,at=x.at)
    axis(2,at=y.at)
    mtext(side=1,line=1.6,las=0,cex=0.7,"Water content (%)")

    # Kleiner et al. (2015) analytical solution from Fig. 4
    lines(k15b$omega,k15b$zeta,col=col[1],lwd=lwd[1],lty=lty[1])

    lines(dat$omega[,kt]*100,dat$zeta,col=col[2],lwd=lwd[2],lty=lty[2])

    legend("topright",bty="n",inset=0.01,col=col,lwd=lwd,lty=lty,c("Kleiner et al. (2015)","Yelmo"))

    box() 
    
    graphics.off()


    #### TIME SERIES ######


    xlim = c(0,1)
    x.at = pretty(xlim,10)
    
    ylim = c(96,108)
    y.at = seq(96,108,by=2)

    myfigure(fldr,"k15expb-ts",asp=1.6,pointsize=12,type=ptype)
    par(xaxs="i",yaxs="i",cex.axis=1.0)

    par(plt=c(0.12,0.95,0.12,0.95))
    
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    abline(h=y.at,v=x.at,col="grey90",lty=1,lwd=0.5)
    axis(1,at=x.at)
    axis(2,at=y.at)
    mtext(side=1,line=1.6,las=0,cex=0.7,"Time (ka)")
    mtext(side=2,line=2.8,las=0,cex=0.7,"Basal enthalpy (kJ kg-1)")

    #abline(h=k15b$enth[1],lwd=0.5,col=col[1],lty=3)
    lines(k15bts$time,k15bts$enth,col=col[1],lwd=lwd[1],lty=lty[1])

    k = which.min(abs(k15b$zeta - 0.05))
    abline(h=k15b$enth[k],lwd=0.5,col=col[1],lty=lty[1])

    kb = 1 
    lines(dat$time,dat$enth[kb,],col=col[2],lwd=lwd[2],lty=lty[2])
    
    k  = which.min(abs(dat$zeta - 0.05))
    lines(dat$time,dat$enth[k,],col=col[2],lwd=lwd[2]*0.5,lty=lty[2])
    text(xlim[2],103.8,pos=2,cex=0.6,"Values for zeta=0.05")

    legend("bottomright",bty="n",inset=0.01,col=col,lwd=lwd,lty=lty,c("Kleiner et al. (2015)","Yelmo"))

    box() 
    
    graphics.off()

}