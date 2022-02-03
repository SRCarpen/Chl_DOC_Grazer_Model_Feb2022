# Fit DOC-Grazer-Chl model:  GitHub example

# (C) Stephen R Carpenter 2 February 2022

rm(list = ls())
graphics.off()

# FUNCTIONS 

# constants
PAR0 = 600 # surface irradiance
eps0 <- 0.0213  # Light extinction by pure water, m-1
epsDOC = 0.0514  # DOC light extinction coef, m2 g-1
epsP <- 0.0177   # Phytoplankton, m2 (mg chl)-1
CChl = 60 # C:Chl mass ratio
CC.chl = 60 # Chl carrying capacity = rmax * shading coef / Competition coef
CC.C = CC.chl*CChl # C carrying capacity
b = 1.e-4 # interference coef for grazing; Scheffer et al. 2008
ghalf = 0.005*CC.C
q = 2 # exponent >= 2

# Function returns negative log likelihood for one-step prediction of chl
# Input: parameters to be estimated, A as chl, Zt m, C is DOC mg/L, N is vector length
OneStep = function(pars,A0vec,A1vec,Ztvec,Cvec,ZBvec,N) {   
  # unpack parameters
  epar = pars #exp(pars)
  rmax = epar[1]
  hmax = epar[2]
  sigma = epar[3]
  
  # simulate
  Ahat = rep(0,N)
  for(i in 1:N) {   # loop over A0 vector
    Chl = A0vec[i]
    A = Chl*CChl  # phyto biomass as C
    ZT = Ztvec[i]
    C = Cvec[i]
    H = ZBvec[i]
    # break the step into 10 steps
    dt = 0.1
    x0 = A  # C units in loop, except for shading
    for(j in 1:10) {
      # light extinction epilimnion
      x.chl = x0/CChl # convert C to Chl for extinction calculation
      eps.pool = epsDOC*C + epsP*x.chl + eps0
      PARmean = (PAR0/ZT)*(1/eps.pool)*(1 - exp(-eps.pool*ZT))
      # phytoplankton dynamics
      cc = rmax*(PARmean/PAR0)/CC.C # competition parameter
      G = ( hmax*H*x0^q/(ghalf^q + x0^q) ) # type 3 grazing for handling & attack setup
      clog = exp(-b*x0)  # interference with grazing at high A
      dx = (rmax*(PARmean/PAR0)*x0 - cc*x0*x0 - G*clog)
      x1 = x0 + (dx*dt)
      x0 = x1
    }
    Ahat[i] = x0/CChl  # save Chl for comparison to data
  }
  devmodel = Ahat-A0vec
  devobs = A1vec - A0vec
  dev = Ahat - A1vec # compare 1 step prediction to observation
  #dev = devmodel - devobs # compare rates
  #NLL = sum(dev*dev) # least squares, sigma estimate will be bad
  NLL = 0.5*N*log(2*pi) + 0.5*N*log(sigma*sigma) + 
    (sum(dev*dev)/(2*sigma*sigma) )
  return(NLL)
}  # end loss function

# Function returns one-step prediction of chlorophyll, A
# Input: parameters to be estimated, A as chl, Zt m, C is DOC mg/L, N is vector length
Yhat = function(pars,A0vec,A1vec,Ztvec,Cvec,ZBvec,N) {   
  # unpack parameters
  epar = pars #exp(pars)
  rmax = epar[1]
  hmax = epar[2]
  sigma = epar[3]
  
  # simulate
  Ahat = rep(0,N)
  for(i in 1:N) {   # loop over A0 vector
    Chl = A0vec[i]
    A = Chl*CChl  # phyto biomass as C
    ZT = Ztvec[i]
    C = Cvec[i]
    H = ZBvec[i]
    # break the step into 10 steps
    dt = 0.1
    x0 = A  # C units in loop, except for shading
    for(j in 1:10) {
      # light extinction epilimnion
      x.chl = x0/CChl # convert C to Chl for extinction calculation
      eps.pool = epsDOC*C + epsP*x.chl + eps0
      PARmean = (PAR0/ZT)*(1/eps.pool)*(1 - exp(-eps.pool*ZT))
      # phytoplankton dynamics
      cc = rmax*(PARmean/PAR0)/CC.C # competition parameter
      G = ( hmax*H*x0^q/(ghalf^q + x0^q) ) # type 3 grazing for handling & attack setup
      clog = exp(-b*x0)  # interference with grazing at high A
      dx = (rmax*(PARmean/PAR0)*x0 - cc*x0*x0 - G*clog)
      x1 = x0 + (dx*dt)
      x0 = x1
    }
    Ahat[i] = x0/CChl  # save Chl for comparison to data
  }
  return(Ahat)
}  # end one-step projection

# END MODEL FUNCTIONS =================================================

# Input data
# variates: DoY Year darkBGA dark.lBGA  dDOsat  dpH 
#   Chl_Manual Zmix_daily      DOC  ZB_gm2  ZB_gm3
#save(Rdaily,Tdaily,file='Daily_PeterTuesday2015_all+ZB.Rdata')
load(file='Daily_PeterTuesday2015_all+ZB.Rdata')
dat0 = Rdaily # select lake
RChl = dat0$Chl_Manual
NR = length(RChl)
RChl0 = RChl[1:(NR-1)]
RChl1 = RChl[2:NR]
Rtvec = dat0[,1]
#
dat1 = Tdaily
TChl = dat1$Chl_Manual
NT = length(TChl)
TChl0 = TChl[1:(NT-1)]
TChl1 = TChl[2:NT]
Ttvec = dat1[,1]
#
# Merge Chl
Chl0 = c(RChl0,TChl0)
Chl1 = c(RChl1,TChl1)
tvec = c(Rtvec,Ttvec)
tplot = c(Rtvec[2:NR],Ttvec[2:NR])
NChl = NR+NT

windows()
par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(Chl0,Chl1,type='p',pch=19,cex=0.7,col='ForestGreen',xlab='A0',ylab='A1')
abline(a=0,b=1,lwd=2,lty=2,col='black')

# Set up for fitting - covariates
RZtvec = dat0$Zmix_daily[1:(NR-1)]
RDOCvec = dat0$DOC[1:(NR-1)]
RZBvec = dat0$ZB_gm3 * 0.5 * 1000 # dry -> C and mg/L -> ug/L
print('Peter mean ZT, DOC, Chl, ZB',quote=F)
print(c(mean(RZtvec),mean(RDOCvec),mean(RChl0),mean(RZBvec)),quote=F)
#
TZtvec = dat1$Zmix_daily[1:(NR-1)]
TDOCvec = dat1$DOC[1:(NR-1)]
TZBvec = dat1$ZB_gm3 * 0.5 * 1000 # dry -> C and mg/L -> ug/L
print('Tuesday mean ZT, DOC, Chl, ZB',quote=F)
print(c(mean(TZtvec),mean(TDOCvec),mean(TChl0),mean(TZBvec)),quote=F)
#
# Merge covariates
Ztvec = c(RZtvec,TZtvec)
DOCvec = c(RDOCvec,TDOCvec)
ZBvec = c(RZBvec,TZBvec)

# Retrieve multi-lake regression model to predict Zt 
load(file='ZmixModel_LakeYears_Hbird.Rdata')  # regression of lake years
print(summary(reg1)) 
pars = reg1$coefficients
#
# Use median squarea sqrt(lake area in ha)
# medsquarea = median(squarea,na.rm=T)  # use this line for median of all data
#medsquarea = median(LYall$Squarea,na.rm=T)  # use this line for median of lake years
medsquarea = sqrt(2.7*1.e4/pi)  # use this line for Peter Lake example


# Plot covariates
windows(width=12,height=12)
par(mfrow=c(4,2),mar=c(2,4.4,2,2)+0.1,cex.axis=1.8,cex.lab=2)
# Chl
plot(Rtvec,RChl,type='l',lwd=2,col='forestgreen',xlab='',
     ylab='Chlorophyll',main='Peter Lake 2015')
plot(Ttvec,TChl,type='l',lwd=2,col='forestgreen',xlab='',
     ylab='Chlorophyll',main='Tuesday Lake 2015')
par(mar=c(3.5,4.5,1,2)+0.1)
# Zoop Biom
plot(Rtvec,RZBvec,type='l',lwd=2,col='blue',xlab='',
     ylab='Zooplankton')
plot(Ttvec,TZBvec,type='l',lwd=2,col='blue',xlab='',
     ylab='Zooplankton')
# DOC
plot(Rtvec[1:(NR-1)],RDOCvec,type='l',lwd=2,col='sienna',xlab='',
     ylab='DOC')
plot(Ttvec[1:(NR-1)],TDOCvec,type='l',lwd=2,col='sienna',xlab='',
     ylab='DOC')
# Thermocline depth
par(mar=c(4,4.5,1,2)+0.1)
plot(Rtvec[1:(NR-1)],RZtvec,type='l',lwd=2,col='black',xlab='Day of Year',
     ylab='Mixed Depth')
plot(Ttvec[1:(NR-1)],TZtvec,type='l',lwd=2,col='black',xlab='Day of Year 2015',
     ylab='Mixed Depth')

# guess parameters rmax, hmax, sigma
sdhat = sd( (Chl1-Chl0),na.rm=T)
guess = c(3,0.2,sdhat)
#guess = log( c(3,0.2,sdhat) ) 
# fit model
fit1 = optim(guess,OneStep,gr=NULL,Chl0,Chl1,Ztvec,DOCvec,ZBvec,(NChl-2),
             method='Nelder-Mead')
             #method='SANN')
print('Optim results',quote=F)
print(c('parameters',fit1$par),quote=F)
#print(c('parameters',exp(fit1$par)),quote=F)
print(c('final NLL',fit1$value),quote=F)
print(c('iterations',fit1$counts),quote=F)
print(c('convergence, 0 is good',fit1$convergence),quote=F)

Ahat = Yhat(fit1$par,Chl0,Chl1,Ztvec,DOCvec,ZBvec,(NChl-2))

windows(width=12,height=4)
par(mfrow=c(1,3),mar=c(4,4.5,2,2)+0.1,cex.axis=1.6,cex.lab=1.8)
yrange=range(Ahat)
xrange=range(Chl1)
plot(Ahat[1:(NR-1)],Chl1[1:(NR-1)],type='p',pch=19,cex=1.2,ylim=yrange,
     xlim=xrange,col='ForestGreen',xlab='Predicted Chl',
     ylab='Observed Chl')
points(Ahat[NR:(NR+NT-2)],Chl1[NR:(NR+NT-2)],type='p',pch=17,cex=1.2,
       col='sienna')
abline(a=0,b=1,lwd=2,lty=2,col='black')
yrange=range(Ahat,Chl1)
plot(tplot[1:(NR-1)],Ahat[1:(NR-1)],type='l',lwd=2,col='blue',ylim=yrange,
     xlab='DoY',ylab='Chlorophyll',main='Peter Lake')
points(tplot[1:(NR-1)],Chl1[1:(NR-1)],type='p',pch=19,cex=1,col='ForestGreen')
#
plot(tplot[NR:(NR+NT-2)],Ahat[NR:(NR+NT-2)],type='l',lwd=2,col='red',ylim=yrange,
     xlab='DoY',ylab='Chlorophyll',main='Tuesday Lake')
points(tplot[NR:(NR+NT-2)],Chl1[NR:(NR+NT-2)],type='p',pch=19,cex=1,col='sienna')

print('correlation of Chl with one-step observation & prediction',quote=F)
print(cor(Chl0,Chl1))
print(cor(Chl1,Ahat))

# Calculate "no model" correlation of A1 and A0
N = length(Chl1)
dev = Chl1-Chl0
sigma = sd(dev)
NLL0 = 0.5*N*log(2*pi) + 0.5*N*log(sigma*sigma) + 
  (sum(dev*dev)/(2*sigma*sigma) )
NLLmod = fit1$value
print('No model NLL and model NLL',quote=F)
print(c(NLL0,NLLmod))
AICmod = 2*3 + 2*NLLmod
AIC0 = 2*2 + 2*NLL0

# Calculate linear model relation of A1 and A0
linA1A0 = lm(Chl1 ~ Chl0)
print('summary of linear model A1 ~ A0',quote=F)
print(summary(linA1A0))
lindev = linA1A0$residuals
lsigma = sd(lindev)
NLL1 = 0.5*N*log(2*pi) + 0.5*N*log(lsigma*lsigma) + 
  (sum(lindev*lindev)/(2*lsigma*lsigma) )
NLLmod = fit1$value
AICmod = 2*3 + 2*NLLmod
AIC1 = 2*3 + 2*NLL1

print('No model, Linear model NLL and model NLL',quote=F)
print(c(NLL0,NLL1,NLLmod))

print('No model AIC, linear model AIC, and model AIC',quote=F)
print(c(AIC0,AIC1,AICmod))

windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4,4.5,2,2)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(Chl0,Chl1,type='p',pch=19,cex=0.7,col='ForestGreen',xlab='A0',ylab='A1')
abline(a=0,b=1,lwd=2,lty=2,col='black')
plot(Chl1-Chl0,Ahat-Chl0,type='p',pch=19,cex=0.7,col='ForestGreen',xlab='Obs Rate',
     ylab='Model Rate')
abline(a=0,b=1,lwd=2,lty=2,col='black')

# ==========================================================================
# Plot functions from model ================================

# Define lake-specific terms; note 
# (mean(RZtvec),mean(RDOCvec),mean(RChl0),mean(RZBvec)
RDOCf = mean(RDOCvec) # Choose DOC for lake to be simulated
RZT = mean(RZtvec)
RZbiom = mean(RZBvec) # Lake mean  #0.08*CC.C # from parameters

TDOCf = mean(TDOCvec) # Choose DOC for lake to be simulated
TZT = mean(TZtvec)
TZbiom = mean(TZBvec) # Lake mean  #0.08*CC.C # from parameters

PAR0 <- 600  # Surface irradiance, microEinsteins m-1 s-1; 600 is typical daily median
# Light extinction parameters, Carpenter et al. L&O 1998
eps0 <- 0.0213  # Light extinction by pure water, m-1
epsDOC = 0.0514  # DOC light extinction coef, m2 g-1
epsP <- 0.0177   # Phytoplankton, m2 (mg chl)-1

# Phytoplankton dynamics

rmax = fit1$par[1] #3 # max unshaded growth coef, nominal 3
CC.chl = 60 # Chl carrying capacity = rmax * shading coef / Competition coef
CC.C = CC.chl*CChl # C carrying capacity

# Grazing function
hmax = fit1$par[2] # fitted attack rate

# Test growth and grazing curves

# general phytoplankton growth & grazing function
# Phytoplankton growth function
grograze = function(C0,p0,ZT,Zb) { 
  peffect = p0/5  # p effect on CC for test case; load / high experimental load
  Agrad = seq(0,CC.C*peffect*1.1,length.out=200)
  Agro = rep(0,200)
  for(i in 1:200) {
    Chl = Agrad[i]/CChl  # convert to chlorophyll to compute shading
    eps.pool = epsDOC*C0 + epsP*Chl + eps0 
    PARmean = (PAR0/ZT)*(1/eps.pool)*(1 - exp(-eps.pool*ZT))
    cc = rmax*(PARmean/PAR0)/(CC.C*peffect) # competition parameter
    Agro[i] = (rmax*(PARmean/PAR0)*Agrad[i] - cc*Agrad[i]^2) 
  }
  # Grazing loss function
  G = ( hmax*Zb*Agrad^q/(ghalf^q + Agrad^q) )*exp(-b*Agrad)
  # Net growth
  Net = Agro - G
  # Equilibria for plotting
  srate = sign(Agro-G)
  dsrate = diff(srate)
  ieq = which(dsrate!=0)
  Aeq.C = Agrad[ieq]
  # List to return
  outlist = list(Agrad,Agro,G,Net,Aeq.C)
  return(outlist)
}

# Calculate growth curves and eq for Peter Lake
p0 = 3 # choose P load
Rout = grograze(RDOCf,p0,RZT,RZbiom)
RAgrad = Rout[[1]]
Rgro = Rout[[2]]
Rgrz = Rout[[3]]
Rnet = Rout[[4]]
Req.C = Rout[[5]]
Req.chl = Req.C/CChl

# Plots for Peter Lake
windows(width=5,height=10)
par(mfrow=c(2,1),mar=c(2.5,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
yrange = range(Rgro,Rgrz)
plot(RAgrad,Rgrz,type='l',lwd=2,col='blue',ylim=yrange,
     xlab='Phytoplankton (C)',ylab='Growth & Grazing',
     main=paste('DOC = ',round(RDOCf,1),', P load = ',p0) )
points(RAgrad,Rgro,type='l',lwd=2,col='forestgreen')
abline(h=0,lty=3,lwd=2)
#
par(mar=c(4,4.5,1,2)+0.1)
plot(RAgrad,Rnet,type='l',lwd=2,col='blue',
     xlab='Phytoplankton (C)',ylab='Growth - Grazing') 
#main=paste('DOC = ',round(DOCf,1),', P load = ',peffect*5) )
abline(h=0,lty=3,lwd=2)
points(Req.C,c(0,0,0),type='p',pch=c(19,21,19),lwd=2,cex=2.5,col='black')

# Calculate growth curves and eq for Tuesday Lake 
p0 = 3 # choose P load
Tout = grograze(TDOCf,p0,TZT,TZbiom)
TAgrad = Tout[[1]]
Tgro = Tout[[2]]
Tgrz = Tout[[3]]
Tnet = Tout[[4]]
Teq.C = Tout[[5]]
Teq.chl = Teq.C/CChl

# Plots for Tuesday Lake
windows(width=5,height=10)
par(mfrow=c(2,1),mar=c(2.5,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
yrange = range(Tgro,Tgrz)
plot(TAgrad,Tgrz,type='l',lwd=2,col='blue',ylim=yrange,
     xlab='Phytoplankton (C)',ylab='Growth & Grazing',
     main=paste('DOC = ',round(TDOCf,1),', P load = ',p0) )
points(TAgrad,Tgro,type='l',lwd=2,col='forestgreen')
abline(h=0,lty=3,lwd=2)
#
par(mar=c(4,4.5,1,2)+0.1)
plot(TAgrad,Tnet,type='l',lwd=2,col='blue',
     xlab='Phytoplankton (C)',ylab='Growth - Grazing') 
main=paste(c('DOC = ',round(TDOCf,1),', P load = ',p0) )
abline(h=0,lty=3,lwd=2)
points(Teq.C,c(0,0,0),type='p',pch=c(19,21,19),lwd=2,cex=2.5,col='black')

# NUMERICAL SEARCH FOR CRITICAL P LOAD OVER A RANGE OF DOC ==================

# Test derivative method for critical point: 
#  dA/dt = f(A), at a critical point f(A) = 0 and f'(A) = 0
dAdt = function(A.C,cvec) {
  C0 = cvec[1]
  p0 = cvec[2]
  ZT = cvec[3]
  Zbiom = cvec[4]
  # growth accounting for DOC & P load
  peffect = p0/5
  Chl = A.C/CChl  # convert to chlorophyll to compute shading
  eps.pool = epsDOC*C0 + epsP*Chl + eps0
  PARmean = (PAR0/ZT)*(1/eps.pool)*(1 - exp(-eps.pool*ZT))
  cc = rmax*(PARmean/PAR0)/(CC.C*peffect) # competition parameter
  Agro = (rmax*(PARmean/PAR0)*A.C - cc*A.C^2)
  # Grazing
  G = ( hmax*Zbiom*A.C^q/(ghalf^q + A.C^q) )*exp(-b*A.C)
  # net and roots
  net = Agro - G
  return(net)
}

# Function for critical criterion
crit = function(A.C,cvec) {
  C0 = cvec[1]
  p0 = cvec[2]
  ZT = cvec[3]
  Zbiom = cvec[4]
  cvec = c(C0,p0,ZT,Zbiom)
  rate = dAdt(A.C,cvec)
  Aup = A.C*1.01
  Adn = A.C*0.99
  rate.up = dAdt(Aup,cvec)
  rate.dn = dAdt(Adn,cvec)
  drate = ( (rate.up-rate)+(rate-rate.dn) )/2
  dA = (Aup-Adn)/2
  dratedA = drate/dA
  # Note numDeriv grad() fails within optim()
  #dratedA = grad(dAdt,x=A.C,method="Richardson",cvec)
  loss = rate^2 + dratedA^2 # loss function
  return(loss)
}

# Find critical point for Peter Lake at a P load
cvec = c(RDOCf,3,RZT,RZbiom) # vector of constants for the analysis
print('',quote=F)
critest = optim(100,crit,gr=NULL,cvec,
                method='Brent',lower=10,upper=1000)
print('test of numerical critical point finder',quote=F)
print(critest$par)
print(critest$convergence)
print('Equilibria by graphical method',quote=F)
print(Req.C)
print('',quote=F)

# Function to find carbon equilibria by graphical method
# Note that polynomial method will not work because of embedded shading function
modeqC = function(C0,p0,ZT,Zbiom) {  
  # growth accounting for DOC & P load
  peffect = p0/5
  Agrad = seq(0,CC.C*2,length.out=1000)
  Chl = Agrad/CChl  # convert to chlorophyll to compute shading
  eps.pool = epsDOC*C0 + epsP*Chl + eps0
  PARmean = (PAR0/ZT)*(1/eps.pool)*(1 - exp(-eps.pool*ZT))
  cc = rmax*(PARmean/PAR0)/(CC.C*peffect) # competition parameter
  Agro = (rmax*(PARmean/PAR0)*Agrad - cc*Agrad^2)
  # Grazing
  G = ( hmax*Zbiom*Agrad^q/(ghalf^q + Agrad^q) )*exp(-b*Agrad)
  # net and roots
  net = Agro - G
  srate = sign(net)
  dsrate = diff(srate)
  ieq = which(dsrate!=0)
  Aeq.C = Agrad[ieq]
  return(Aeq.C)
}
#

# Study equilibrium on P gradient for low and high DOC ====================
DOCvec = c(1:20)
NC = length(DOCvec)
NP = 100 # number of P loads
pvec = seq(0.3,5,length.out=NP)
neqvec = rep(0,NP)  # disposable vector of number of equilibria
Pcrit = rep(0,NC) # vector to hold P values at flip
#
for(i in 1:NC) {
  for(j in 1:NP) {
    C0 = DOCvec[i]
    p0 = pvec[j]
    ZThat = pars[1] + pars[2]*medsquarea + pars[3]*C0 + pars[4]*p0
    ZThat = max(c(ZThat,0.1))
    eqest = modeqC(C0,p0,ZThat,RZbiom)
    neq = length(eqest)
    neqvec[j] = neq
  }
  i3 = which(neqvec==3)
  Pcrit[i] = min(pvec[i3])
}

windows()
par(mfrow=c(1,1),mar=c(4,4.5,2,1)+0.1,cex.axis=2,cex.lab=2,cex.main=1.5)
plot(DOCvec,Pcrit,type='b',lwd=2,pch=19,cex=1,xlab='DOC',ylab='Critical value')

# Composite 4 panel figure for paper
windows(width=16,height=4)
par(mfrow=c(1,4),mar=c(4,4.5,2,1)+0.1,cex.axis=2,cex.lab=2,cex.main=1.5)
yrange=range(Ahat,Chl1)
plot(tplot[1:(NR-1)],Ahat[1:(NR-1)],type='l',lwd=2,col='blue',ylim=yrange,
     xlab='Day',ylab='Chlorophyll',main='A. Peter Lake')
points(tplot[1:(NR-1)],Chl1[1:(NR-1)],type='p',pch=19,cex=1,col='ForestGreen')
legend('topleft',legend=c('data','model'),pch=c(19,NA),lwd=c(NA,2),cex=1.5,
       col=c('Forestgreen','blue'),bty='n')
#
plot(tplot[NR:(NR+NT-2)],Ahat[NR:(NR+NT-2)],type='l',lwd=2,col='red',ylim=yrange,
     xlab='Day',ylab='Chlorophyll',main='B. Tuesday Lake')
points(tplot[NR:(NR+NT-2)],Chl1[NR:(NR+NT-2)],type='p',pch=19,cex=1,col='sienna')
legend('topleft',legend=c('data','model'),pch=c(19,NA),lwd=c(NA,2),cex=1.5,
       col=c('sienna','red'),bty='n')
#
plot(RAgrad/CChl,Rnet,type='l',lwd=2,col='blue',
     xlab='Chlorophyll, ug/L',ylab='Growth - Grazing',
     main=paste('C.  DOC = ',round(RDOCf,1),', P load = ',3) )
abline(h=0,lty=3,lwd=2)
points(Req.chl,c(0,0,0),type='p',pch=c(19,21,19),lwd=2,cex=2.5,col='black')
#
plot(TAgrad/CChl,Tnet,type='l',lwd=2,col='blue',
     xlab='Chlorophyll, ug/L',ylab='Growth - Grazing',
     main=paste('D.  DOC = ',round(TDOCf,1),', P load = ',3) )
abline(h=0,lty=3,lwd=2)
points(Teq.chl,c(0,0,0),type='p',pch=c(19,21,19),lwd=2,cex=2.5,col='black')

