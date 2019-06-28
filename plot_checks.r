library(myr)

rho_ice  = 910.0 
rho_w    = 1000.0 
cp       = 2009.0           # [J kg-1 K-1]
L_ice    = 333500.0         # [J kg-1] Latent heat

# For comparison with Kleiner et al. (2015), use this T_ref offset for enthalpy==0.0
T_ref = 223.15              # [K] 
E_ref = (T_ref*cp)*1e-3     # [kJ kg-1]

# Set this flag to True if data from icetemp will be in units of Celcius to start with     
is_celcius = FALSE 
    

load_icetemp = function(filename)
{
    dat = my.read.nc(filename)
    dat$filename = filename 
    dat$time = dat$time*1e-3     # [a] => [ka]
    #dat$dz   = dat$H_ice[1]/(length(dat$zeta)-1)                # [m] vertical grid resolution 
    dat$dz   = round(mean(diff(dat$zeta_ac))*dat$H_ice[1],2)    # [m] vertical grid resolution, nominal if unevenly spaced
    
    # Calculate melt rate [m/a w.e.]
    dat$bmb_we = dat$bmb*rho_ice/rho_w 

    # Adjust enthalpy for reference
    dat$enth = dat$enth*1e-3 - E_ref 

    if (! is_celcius) {
        # Convert to Celcius for plots 
        T0 = 273.15 
        dat$T_ice   = dat$T_ice   - T0 
        dat$T_pmp   = dat$T_pmp   - T0 
        dat$T_srf   = dat$T_srf   - T0 
        dat$T_robin = dat$T_robin - T0     
    }

    return(dat)
}

# Load data 
if (TRUE) {

    experiment = "k15expb"

    # Load data from Exp. A
    dat_k15expa = load_icetemp(paste0("test_k15expa.nc"))

    k15expb_filenames = c(
    "test_k15expb_cr0.10E+00_dz0.10E+02.nc",
    "test_k15expb_cr0.10E+00_dz0.50E+00.nc",
    "test_k15expb_cr0.10E-01_dz0.10E+02.nc",
    "test_k15expb_cr0.10E-01_dz0.50E+00.nc",
    "test_k15expb_cr0.10E-02_dz0.10E+02.nc",   # qref2
    "test_k15expb_cr0.10E-02_dz0.50E+00.nc",   # qref1
    "test_k15expb_cr0.10E-03_dz0.10E+02.nc",
    "test_k15expb_cr0.10E-03_dz0.50E+00.nc") 

    n_expb = length(k15expb_filenames)

    dats_k15expb = list() 
    for (q in 1:n_expb) {
        dats_k15expb[[q]] = load_icetemp(k15expb_filenames[q])
    }

    qref1 = 6
    qref2 = 5 

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

    dat = dat_k15expa

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

    dat  = dats_k15expb[[qref1]] 
    dat2 = dats_k15expb[[qref2]] 

    kt  = length(dat$time) 
    kt2 = length(dat2$time) 

    ylim = c(0,1)
    y.at = seq(0,1,by=0.1)

    col = c("#ef8a62","black","#99d8c9")          # reddish color: "#d6604d"
    lwd = c(4.0,2.5,1.6)
    lty = c(1,1,1)

    myfigure(fldr,"k15expb",asp=2.3,pointsize=12,type=ptype)
    par(xaxs="i",yaxs="i",cex.axis=1.0)

    par(mfrow=c(1,3))

    # Panel 1: Enthalpy #############################
    par(plt=c(0.12,0.99,0.12,0.95))
    xlim = c(92,108)
    x.at = c(92,96,100,104,108)
    
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    abline(h=y.at,v=x.at,col="grey90",lty=1,lwd=0.5)
    axis(1,at=x.at)
    axis(2,at=y.at)
    mtext(side=1,line=1.6,las=0,cex=0.7,"Enthalpy (kJ kg-1)")

    # Kleiner et al. (2015) analytical solution from Fig. 4
    lines(k15b$enth,k15b$zeta,col=col[1],lwd=lwd[1],lty=lty[1])

    lines(dat$enth[,kt],dat$zeta,  col=col[2],lwd=lwd[2],lty=lty[2])
    lines(dat2$enth[,kt2],dat2$zeta,col=col[3],lwd=lwd[3],lty=lty[3])

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

    lines(dat$T_ice[,kt],dat$zeta,  col=col[2],lwd=lwd[2],lty=lty[2])
    lines(dat2$T_ice[,kt2],dat2$zeta,col=col[3],lwd=lwd[3],lty=lty[3])
    
    box() 
    
    # Panel 3: omega #############################
    par(plt=c(0.08,0.95,0.12,0.95))
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
    lines(dat2$omega[,kt2]*100,dat2$zeta,  col=col[3],lwd=lwd[3],lty=lty[3])
    
    legend("topright",bty="n",inset=0.01,col=col,lwd=lwd,lty=lty,c("Kleiner et al. (2015)","Yelmo dz=0.5m","Yelmo dz=10m"))

    box() 
    
    ## Panel 3 inset: H_cts ######
    par(new=TRUE,plt=c(0.45,0.90,0.40,0.75),cex.axis=0.8)
    xlim = c(0.8e-4,2e-1)
    x.at = c(1e-4,1e-3,1e-2,1e-1)
    ylim = c(15,42)
    y.at = c(20,25,30,35,40)

    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE,log="x")
    rect(xlim[1],ylim[1],xlim[2],ylim[2],col="white",bg="white")
    #abline(h=y.at,v=x.at,col="grey90",lty=1,lwd=0.5)
    axis(1,at=x.at)
    axis(2,at=y.at)
    mtext(side=1,line=1.2,las=0,cex=0.5,"K0/Kc")
    mtext(side=2,line=1.3,las=0,cex=0.5,"CTS height (m)")
    
    # Analytical solution (roughly)
    abline(h=18.9,col=col[1],lwd=lwd[1],lty=lty[1])

    for (q in 1:length(dats_k15expb)) {
        tmp = dats_k15expb[[q]]
        kt_now = length(tmp$time)
        col_now = "grey80"
        if (tmp$dz==0.5) col_now = col[2]
        if (tmp$dz > 8)  col_now = col[3]
        points(tmp$enth_cr,tmp$H_cts[kt_now],col=col_now,pch=20,cex=1.5)
        #cat(tmp$enth_cr,"  ",tmp$H_cts[kt_now],"\n")
    }
    
    points(dat$enth_cr,dat$H_cts[kt],col=col[2],pch=20,cex=1.4)
    points(dat2$enth_cr,dat2$H_cts[kt2],col=col[3],pch=20,cex=1.4)
    
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
    text(xlim[2],103.6,pos=2,cex=0.6,"Values for zeta=0.05")

    legend("bottomright",bty="n",inset=0.01,col=col,lwd=lwd,lty=lty,c("Kleiner et al. (2015)","Yelmo"))

    box() 
    
    graphics.off()

}