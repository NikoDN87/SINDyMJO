#########################################################################################
#                         Code for applying SINDY to MJO                                #
#                                NDN 08/09/23                                           #
#---------------------------------------------------------------------------------------#
#                                                                                       #
# Idea: 1. Take random samples of MJO's OMI data (first two principal components)  and  #
#          look for bi-dimensional dynamical systems over a polinomial dictionary       #
#       2. Select models using Akaike Information Criteria                              #
#       3. Inferre models according to El Niño and La Niña years                        #
#########################################################################################

#-----------------------------------#
#------- Required Library ----------#
#-----------------------------------#----------------------------------------------------
library(ggplot2)
library(ggdendro)
library(reshape2)
library(ggpubr)
library(R1magic)
library(svDialogs)
library(latex2exp)
library(reshape2)
library(viridis)
#--------------------------#
#--- Working  Directory ---#
#--------------------------#-------------------------------------------------------------
wd = "/home/nikolaos/Documentos/GitHub"
setwd(wd)
source("MJO_polar.R")       # call MJO.polar function
#------------------------#
# Select RMM or OMI data:#
#------------------------#---------------------------------------------------------------
index = dlgInput("Select MJO Index: 'RMM', 'OMI' ")$res
if(index == "RMM") {
  MJO_raw = read.csv("RMM_MJO_INDEX_1974_2022.csv", header=T)
} else if (index == "OMI") {
  MJO_raw = read.csv("OMI_MJO_INDEX_1979_2022.csv", header=T)
  colnames(MJO_raw) = c("Año","Mes","Dia","RMM1","RMM2","Amplitud")
}
#--------------#
# Prepare data #
#--------------#-------------------------------------------------------------------------
init_row = 1                            # day 1 of data
fin_row = which(MJO_raw$Año == 2021 & MJO_raw$Mes == 4 & MJO_raw$Dia == 4)  # last day of data

if(index == "OMI") {
  PC1 = MJO_raw$RMM1    # select OMI PC1
  PC2 = MJO_raw$RMM2    # select OMI PC2
  
  MJO_raw$RMM1 = PC2    # change OMI's PCs to coincide with RMM Pcs
  MJO_raw$RMM2 = -PC1   
  
  rm(PC2, PC1)          # remove PC1 and PC2
}

MJO_raw2 = MJO_raw[seq(init_row, fin_row),]             # select MJO data between day 1 and last day
MJO = MJO_raw2[!(MJO_raw2$Mes==2 & MJO_raw2$Dia==29) ,] # remove february the 29th
añoi = MJO$Año[1];                                      # first year
añof = MJO$Año[length(MJO$Año)]                         # last year

all = seq(añoi,añof)    # all years within period
neutros = c(1980,1981,1985,1989,1990,1992,1993,1996,2001,2003,2012,2013,2019) # ENSO neutral years
niñas = c(1983,1984,1988,1995,1998,1999,2000,2005,2007,2008,2010,2011,2016,2017,2020) # ENSO Niña years
niños = c(1979,1982,1986,1987,1991,1994,1997,2002,2004,2006,2009,2014,2015,2018) # ENSO Niño years

# select years according to ENSO phases or full period
años.ch = dlgInput("Select years: 'all', 'neutros', 'niñas', 'niños' ")$res
if(años.ch == "all") {años = all} else if (años.ch == "neutros") {años = neutros
} else if (años.ch == "niñas") {años = niñas} else if (años.ch == "niños") {años = niños}

x_in.df = data.frame(vector())           # empty data.frame
for(año in años[1:(length(años))]) {
  
  if(año!=2021) {
    if(index == "RMM") {
      # Select days corresponding to DJFM months (and last [first] days of november [april])
      x_in.0 = MJO[(MJO$Año==año & MJO$Mes==11 & MJO$Dia==27) |     # las days of november
                     (MJO$Año==año & MJO$Mes==11 & MJO$Dia==28) |
                     (MJO$Año==año & MJO$Mes==11 & MJO$Dia==29) |
                     (MJO$Año==año & MJO$Mes==11 & MJO$Dia==30) |
                     (MJO$Año==año & MJO$Mes==12) |                 # december
                     ((MJO$Año==(año+1) & MJO$Mes==1)) |            # january
                     (MJO$Año==(año+1) & MJO$Mes==2) |              # february
                     (MJO$Año==(año+1) & MJO$Mes==3) |              # march
                     (MJO$Año==(año+1) & MJO$Mes==4 & MJO$Dia==1) | # first days of april
                     (MJO$Año==(año+1) & MJO$Mes==4 & MJO$Dia==2) |
                     (MJO$Año==(año+1) & MJO$Mes==4 & MJO$Dia==3) |
                     (MJO$Año==(año+1) & MJO$Mes==4 & MJO$Dia==4),]
      
      # Apply 9-day moving average to smooth MJO's RMM index
      x_in = x_in.0[5:125,]                               # select days indexes corresponding to DJFM months
      for(t in 1:length(x_in[,1])) {
        x_in$RMM1[t] =  mean(x_in.0$RMM1[seq(t,t+8)])  # 9-day moving average centerd in day 5
        x_in$RMM2[t] =  mean(x_in.0$RMM2[seq(t,t+8)])
      }
    } else if (index == "OMI") { # OMI does not requiere smoothing
      x_in = MJO[(MJO$Año==año & MJO$Mes==12) |                   # december
                   ((MJO$Año==(año+1) & MJO$Mes==1)) |            # january
                   ((MJO$Año==(año+1) & MJO$Mes==2)) |            # february
                   ((MJO$Año==(año+1) & MJO$Mes==3)),]            # march
    }
    chg.year = which(x_in$Año==año & x_in$Mes==12)
    x_in$Año[chg.year] = x_in$Año[chg.year]+1        # change year corresponding to december, 
                                                     # in order to coincide with year of following months
    ll = length(x_in$RMM1)
    vel1 = rep(NA, length(x_in$RMM1))    # empty vector for velocity component 1
    vel2 = rep(NA, length(x_in$RMM1))    # empty vector for veloity component 2
    # numericaly calculate velocity with 3 and 4 points schemes
    for(t in 1:length(x_in[,1])) {
      if ( t <=3) {
        vel1[t] = (-3*x_in$RMM1[t] + 4*x_in$RMM1[t+1] - x_in$RMM1[t+2])/2
        vel2[t] = (-3*x_in$RMM2[t] + 4*x_in$RMM2[t+1] - x_in$RMM2[t+2])/2
      }
      else if (t>=(length(x_in$RMM1)-2)) {
        vel1[t] = (3*x_in$RMM1[t] - 4*x_in$RMM1[t-1] + x_in$RMM1[t-2])/2
        vel2[t] = (3*x_in$RMM2[t] - 4*x_in$RMM2[t-1] + x_in$RMM2[t-2])/2
      }
      else {
        vel1[t] = (-x_in$RMM1[t+2] + 8*x_in$RMM1[t+1] - 8*x_in$RMM1[t-1] + x_in$RMM1[t-2])/12
        vel2[t] = (-x_in$RMM2[t+2] + 8*x_in$RMM2[t+1] - 8*x_in$RMM2[t-1] + x_in$RMM2[t-2])/12
      }
    }
    
    # add rows with position and velocity data.frame for year año, to previous x_in.df
    x_in.df = rbind(x_in.df, 
                    cbind(x_in$Año, x_in$RMM1,x_in$RMM2, vel1,vel2))
  }
}
colnames(x_in.df) = c("Año","RMM1","RMM2","RMM1.vel","RMM2.vel")  # column names
# If index = RMM , we need to remove year 1978-1979 (missing data).
if(index == "RMM") {
  rm.index = which(x_in.df$Año==1979 | x_in.df$Año==1978)
  if(length(rm.index)>=1) {x_in.df = x_in.df[-rm.index,]} # RMM no tiene todos los datos del año 78-79
}

# Use MJO.polar to calculate MJO's daily phase and amplitude
DEFM.polar = MJO.polar(x_in.df$RMM1, x_in.df$RMM2)
DEFM.df = cbind(x_in.df, DEFM.polar)              # add columns with amplitude and phase
Amp.ind = which(DEFM.df$Amplitude>1)              # rows with amplutde>1
DEFM.amp.df = DEFM.df[Amp.ind,]                   # select data with amplitude>1

#---------------------------------------#
# Perform random sampling of MJO's data #
#---------------------------------------#------------------------------------------------
MJO.tray = vector("list")                         # empty list
i=0                                               # initiate index i = 0
Nreal = 1000                                      # number of realisations
Nsamp = 2^8                                       # number of random samples/realisation
for(chi in 1:Nreal) {
  
  i=i+1
  # Randomly sample -- without replacement -- over PC1-PC2 phase space
  samp.ind = sample(x = seq(1,length(DEFM.amp.df$RMM1)), size = Nsamp, replace = F  )
  MJO.tray[[i]] = DEFM.amp.df[samp.ind,]  # save sampled data for realization i
  
}
names(MJO.tray) = paste0("Samp",seq(1,Nreal))    # name each data.frame realisation
#-------------#
# Apply SINDy #
#-------------#--------------------------------------------------------------------------
Sindy = function(Psi, xpunto, init.guess, L1.lambda){
  Psi.dim = dim(Psi)
  if(is.null(Psi.dim)) {T.mat = 1} else {T.mat = diag(length(Psi[1,]))}
  out = solveL1(Psi, y = xpunto, T = T.mat, x0 = init.guess , lambda = L1.lambda)$estimate
  return(out)}

# Define parameters and arrays needed for SINDy: #
n = 2                                                   # maximum polynomial power
D = 2                                                   # system dimension
NPsi.0 = (1+n)^D                                        # number of predictors 
NPsi = NPsi.0
Niter = NPsi                                            # number of iterations
Lambda.inv = seq(1,100,1)                               # inverse of hard thresholding parameter lambda
# We select up to 4 models per realisation:
Modelo1.list = vector("list", length = length(MJO.tray)) # inititate empty list for model 1
Modelo2.list = vector("list", length = length(MJO.tray)) # idem model 2
Modelo3.list = vector("list", length = length(MJO.tray)) # idem model 3
Modelo4.list = vector("list", length = length(MJO.tray)) # idem model 4
AIC.list = vector("list", length = length(MJO.tray))     # initiate list for AIC plot/realsation
for(M in 1:length(MJO.tray)) { # M counts realisation
  
  x1.0 = MJO.tray[[M]]$RMM1; x2.0 = MJO.tray[[M]]$RMM2   # position data
  xpunto.0 = cbind(MJO.tray[[M]]$RMM1.vel, MJO.tray[[M]]$RMM2.vel)  # velocity data
  xpunto = xpunto.0
  # Construct matriz Psi with polynomial columns #
  Psi.0 = array(0, dim = c(length(xpunto[,1]),NPsi.0))      # Psi[1,] = [1, xn, yn, xn*yn,.... ]
  for (t in 1:(length(xpunto[,1]))) {                       
    j=0
    for (l1 in 0:n) {
      for(l2 in 0:n) {
        j = j+1
        Psi.0[t,j] = (x1.0[t]^l1)*(x2.0[t]^l2)        # define column corresponding to x1^l1*x2^l2
      }
    }
  }
  # Normalize columns of Psi according to L2-norm #
  L2Psi.0 = rep(0, length(Psi.0[1,]))         # empty vector for L2-norm
  for (j in 1:length(Psi.0[1,])) {
    L2Psi.0[j] = norm(Psi.0[,j],"2")          # j-esima column L2-norm
    Psi.0[,j] = Psi.0[,j]/L2Psi.0[j]          # divide j colum by its L2-norm
  }
  
  Psi = Psi.0     # save normalized matrix Psi.0 in Psi
  L2Psi = L2Psi.0 # idem for vector of L2-norms
  coefs = array(0, dim = c(NPsi,length(Lambda.inv),D))    # initiate null array of models coefficients
  rmse.df = data.frame(Lambda.inv, rep(NA,length(Lambda.inv)),rep(NA,length(Lambda.inv))) # initiate data.frame for calculate rmse for each lambda
  colnames(rmse.df) = c("lambda.inv", "train", "val")
  AIC.df = data.frame(Lambda.inv, rep(NA, length(Lambda.inv)))  # initiate data.frame dor AIC
  colnames(AIC.df) = c("lambda.inv","train")
  for(k in 1:length(Lambda.inv)) {
    for (dim in 1:D) {
      p = qr.solve(Psi, xpunto[,dim])   # initial guess for linear regression problem
      AuxCoefs = Sindy(Psi,xpunto[,dim], init.guess = p, L1.lambda = 0) # apply SINDy 

      for (it in 1:Niter) {
        
        IZ = which(abs(AuxCoefs/L2Psi) < 1/Lambda.inv[k])         # look for coefs < Lambda
        INZ = which(abs(AuxCoefs/L2Psi)>=1/Lambda.inv[k])         # look for coefs >= Lambda
        AuxCoefs[IZ] = 0                                 # impose IZ coefs to be null
        # re-calculate non-null coefficients with reduced linear regression
        if(length(INZ>0)) {
          AuxCoefs[INZ] = Sindy(Psi =  Psi[,INZ], xpunto =  xpunto[,dim], init.guess = AuxCoefs[INZ], L1.lambda = 0 )
        }
        if (it == Niter) {
          coefs[,k,dim] = AuxCoefs/L2Psi} # divide AuxCoefs by L2-norm in order to obtain non-normalise coeffs
      }
    }
    rm(p) # remove initial guess
    #---------------------------------------------#
    # Reconstruct data according to infer model   #--------------------------------------#
    # For each sampled data (x,y), we calculate (x_dot, y_dot) according to the inferred #
    # model, and compare with actual velocity data through rmse and AIC                  #
    # --------------------------------------------#--------------------------------------#
    NCNNx = length(which(coefs[,k,1]!=0))   # number of non-null coefficients for x_dot eq.
    NCNNy = length(which(coefs[,k,2]!=0))   # number of non-null coefficients for y_dot eq.
    NCNN = NCNNx + NCNNy                    # total number of non-null coefficients
    xpunto1.rec = Psi%*%coefs[,k,1]       # componente 1 de modelo
    xpunto2.rec = Psi%*%coefs[,k,2]       # componente 2 de modelo
    
    if(!any(is.nan(xpunto1.rec), is.infinite(xpunto1.rec)) # remove models with possible nan values
       & !any(is.nan(xpunto2.rec), is.infinite(xpunto2.rec)) ) {
      # calculate bi-rmse index
      rmse.df$train[k] = sqrt(sum( (xpunto[,1]-xpunto1.rec)^2 + (xpunto[,2]-xpunto2.rec)^2, 
                                   na.rm = T  )/length(xpunto[,1]))
      # calculate AIC value (least square expression for finite sample)
      AIC.df$train[k] = length(xpunto[,1])*log(rmse.df$train[k]^2) + 2*(NCNN+1) +
        2*(NCNN+1)*(NCNN+2)/(length(xpunto[,1])-NCNN-2)
    }
  }
  # Define DeltaAIC (shift values by minimum)
  AIC.df$train = AIC.df$train - min(AIC.df$train, na.rm = T)  # remove minimum value of AIC
  
  if(sum(!is.na(AIC.df$train))>0) {
    gg = ggplot() +
      geom_rect(aes(ymin = 0, ymax = 2, xmin = -Inf, xmax = Inf), fill = "green3", alpha = 0.3) +
      geom_rect(aes(ymin = 3, ymax = 7, xmin = -Inf, xmax = Inf), fill = "darkorange", alpha = 0.3) +
      geom_point(data = AIC.df, aes(x = lambda.inv, y = train, shape = "sh1", col = "c1")
                 , size = 2.5 ) +
      scale_y_continuous(TeX("$\\Delta$ AIC"), limits = c(0,10)) +
      scale_x_continuous(TeX("$\\lambda^{-1}$"), limits = c(0,100)) +
      theme(axis.title = element_text(size = 16, face = "bold")) +
      theme(axis.text = element_text(size = 15)) + theme_bw()+
      scale_colour_manual(name = "",  values = c("c1"="darkorchid3"),
                          labels = c("Delta AIC")) +
      scale_shape_manual(name = "", values = c("sh1"= 8),
                         labels = c("Delta AIC")) + 
      theme(legend.position = "")

    AIC.list[[M]] = print(gg) # comment line to avoid ploting
  }
  # Select Models #
  AIC.df$train = round(AIC.df$train*10000)/10000
  km1 = which(AIC.df$train == min(AIC.df$train,na.rm = T))     # indexes for AIC min value
  Modelo1 = km1[1]                                             # select first index
  NCNNx1 = length(which(coefs[,Modelo1,1]!=0))                 # Number of non-null coefficients for dot_x   
  NCNNy1 = length(which(coefs[,Modelo1,2]!=0))                 # Idem for dot_y
  km2 = which(AIC.df$train == min(AIC.df$train[-km1], na.rm = T)) # idem for second min value
  Modelo2 = km2[1]                                             # select first index
  NCNNx2 = length(which(coefs[,Modelo2,1]!=0))                 # Number of non-null coefficients for dot_x   
  NCNNy2 = length(which(coefs[,Modelo2,2]!=0))                 # Idem for dot_y
  km3 = which(AIC.df$train == min(AIC.df$train[-c(km1,km2)], na.rm = T)) # idem for third min value
  Modelo3 = km3[1]                                             # select first index
  NCNNx3 = length(which(coefs[,Modelo3,1]!=0))                 # Number of non-null coefficients for dot_x   
  NCNNy3 = length(which(coefs[,Modelo3,2]!=0))                 # Idem for dot_y
  km4 = which(AIC.df$train == min(AIC.df$train[-c(km1,km2,km3)], na.rm = T)) # idem for fourth min value
  Modelo4 = km4[1]                                             # select first index
  NCNNx4 = length(which(coefs[,Modelo4,1]!=0))                 # Number of non-null coefficients for dot_x   
  NCNNy4 = length(which(coefs[,Modelo4,2]!=0))                 # Idem for dot_y
  if(NCNNx1!=0 & NCNNy1!=0) { Modelo1.list[[M]] = coefs[,Modelo1,]  } # select first model
  if(!is.na(Modelo2) &  AIC.df$train[Modelo2]<=5 & NCNNx2!=0 & NCNNy2!=0) {Modelo2.list[[M]] = coefs[,Modelo2,]} # idem model2
  if(!is.na(Modelo3) &  AIC.df$train[Modelo3]<=7 & NCNNx3!=0 & NCNNy3!=0) {Modelo3.list[[M]] = coefs[,Modelo3,]} # idem model 3
  if(!is.na(Modelo4) &  AIC.df$train[Modelo4]<=7 & NCNNx4!=0 & NCNNy4!=0) {Modelo4.list[[M]] = coefs[,Modelo4,]} # idem model 4
}

# Define data.framde with all inferred models:
modelos0.df = data.frame(row.names = seq(1,2*NPsi))
for(i in 1:4) {
  if (i == 1) {  Modelo.list = Modelo1.list} else if (i == 2) {  Modelo.list = Modelo2.list
  } else if (i == 3) {  Modelo.list = Modelo3.list} else if (i == 4) {  Modelo.list = Modelo4.list}
  for(chi in 1:Nreal) {
    if(!is.null(Modelo.list[[chi]])) {
      modelos0.df= cbind(modelos0.df, as.vector(Modelo.list[[chi]]))
    }
  }
}
#------------------------------------#
# Figure: statistics of coefficients:#
#------------------------------------#---------------------------------------------------
num.coefs = rep(NA,2*NPsi)
avg.coefs = rep(NA,2*NPsi)
p75.coefs = rep(NA,2*NPsi)
p25.coefs = rep(NA,2*NPsi)
min.coefs = rep(NA,2*NPsi)
max.coefs = rep(NA,2*NPsi)
for(i in 1:(2*NPsi)) {
  num.coefs[i] = length(which(modelos0.df[i,]!=0))/length(modelos0.df[1,])
  avg.coefs[i] = mean(as.numeric(modelos0.df[i,]), na.rm = T)
  p75.coefs[i] = quantile(as.numeric(modelos0.df[i,]), probs = 0.75)[[1]]
  p25.coefs[i] = quantile(as.numeric(modelos0.df[i,]), probs = 0.25)[[1]]
  min.coefs[i] = min(as.numeric(modelos0.df[i,]), na.rm = T)
  max.coefs[i] = max(as.numeric(modelos0.df[i,]), na.rm = T)
}

avg.modelo = cbind(avg.coefs[1:NPsi], avg.coefs[(NPsi+1):(2*NPsi)])
Coefs.df = data.frame(seq(1,2*NPsi),num.coefs,avg.coefs, p25.coefs, p75.coefs, min.coefs, max.coefs)
colnames(Coefs.df) = c("Coef","Freq","Mean","p25","p75","min","max")

x.lab = rep(TeX( input = c("$1$", "$\\y$", "$\\y^2$", "$\\x$", "$\\xy$", "$\\xy^2$", "$\\x^2$", "$\\x^2y$", "$\\x^2y^2$" )),2)
gg1 = ggplot() + geom_col(data = Coefs.df, aes(x=Coef, y = Freq, fill = ifelse(Coef<=NPsi,"fc1","fc2")),
                          col = "gray35") +
  scale_fill_manual(values = c("fc1"="cornflowerblue", "fc2"="darkolivegreen3"), 
                    labels = c(TeX("$\\dot{x}$"),TeX("$\\dot{y}$")), name = "Equation" ) +
  scale_x_continuous( breaks = seq(1,18), labels = x.lab ) + 
  theme(legend.position = "top", legend.direction = "horizontal") + 
  ylab("Frequency") +
  theme_bw() +
  theme(axis.text.x  = element_text(size = 15, face = "bold", vjust = 0.2),
        axis.title.x = element_text(size = 16, face="bold")) +
  theme(axis.text.y  = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size=16, face="bold")) +
  theme(legend.text = element_text(size = 15, face = "bold"), legend.direction = "horizontal",
        legend.position = "top",
        legend.title = element_text(size = 15, face = "bold")) 

gg2 = ggplot() + geom_col(data = Coefs.df, aes(x=Coef, y = Mean, fill = ifelse(Coef<=NPsi,"fc1","fc2")),
                          col = "gray35") +
  scale_fill_manual(values = c("fc1"="cornflowerblue", "fc2"="darkolivegreen3"), 
                    labels = c(TeX("$\\dot{x}$"),TeX("$\\dot{y}$")), name = "Equation" ) +
  geom_errorbar(data = Coefs.df, aes(x=Coef, ymin=p25, ymax = p75))+
  scale_x_continuous( breaks = seq(1,18), labels = x.lab ) + 
  geom_point(data = Coefs.df, aes(x = Coef, y =min), 
             shape = 25, col = "red3", fill = "white", size = 2) + 
  geom_point(data = Coefs.df, aes(x = Coef, y =max), 
             shape = 24, col = "red3", fill = "white", size = 2) +
  theme_bw() + theme(legend.position = "none") +
  theme(axis.text.x  = element_text(size = 15, face = "bold", vjust = 0.2),
        axis.title.x = element_text(size = 16, face="bold")) +
  theme(axis.text.y  = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size=16, face="bold")) 


ggarrange(nrow = 2, gg1 + rremove("xlab") + rremove("x.text"), gg2)
# to save plot, uncomment the following lines:
# ggsave("Coefs.png", device = png(), height = 6, width = 9,
#       dpi = 300, bg = "white")

#----------------------------------------------------------------------#
# Select different models according to non-null coefficients structure:#
#----------------------------------------------------------------------#-----------------
M.aux1 = modelos0.df[,1]             # select first model among all inferred models
M.aux1[M.aux1!=0] = 1                # identified non-null coefficients (mask for model 1)
mask.df = data.frame(M.aux1)         # initiate mask.df with mask of first model
modelos.df = data.frame(M.aux1)      # initiate modelos.df with coefficients of first model  
for(j in 2:dim(modelos0.df)[2] ) {
  
  aux = modelos0.df[,j]              # select j-ésimo model
  aux[aux!=0] = 1                    # mask for j-ésimo model
  mask.df = cbind(mask.df, aux)      # add in columns, succesive masks 
  L.modelos = dim(modelos.df)[2]     # number of models selected up to iteration j 
  for(l in 1:L.modelos) {
    
    aux0 = modelos.df[,l]            # select mask for model l 
    sum0 = sum(aux0==aux)            # check if mask for model l is equal to mask for model j
    if(sum0 == (2*NPsi)) {break} else {
      if(l<L.modelos) {next} else if(l == L.modelos) { modelos.df = cbind(modelos.df, aux) }   # guardo únicamente si es distinto a todos los ya guardados
    }
  }
}
L.modelos = dim(modelos.df)[2]   # number of models with different structures

# Calculte frequency of occurrence of each of the different models found
avg.list = vector("list")
freq.modelos = rep(NA, L.modelos)
linear.entries = c(1,2,4,10,11,13)          # linear entries of the system
non.linear.entries = seq(1,2*NPsi)[-linear.entries]   # non-linear entries
omega2 = numeric(0)                         # initiate empty value for square angular frequency
for(M in 1:L.modelos) { # M counts over different models found
  
  mask.aux = modelos.df[,M]        # mask for M model
  out = apply(mask.df, MARGIN = 2, function(x)  ifelse( (sum(x==mask.aux) == 2*NPsi),1,0) )
  if(length(which(out==1)) == 1) { #length(which(out==1)) counts number of models with mas.aux structure
    out.avg = modelos0.df[,out==1]
  } else { 
    out.avg = apply(modelos0.df[,which(out==1)], MARGIN = 1, mean)
  }
  
  avg.list[[M]] = cbind(out.avg[1:NPsi], out.avg[(NPsi+1):(2*NPsi)])
  freq.modelos[M] = length( which(out==1) )/length(modelos0.df)
  
  # If all non-linear entries of mas.aux are null, calculate omega2
  if(!any(mask.aux[non.linear.entries]!=0)) {
    if(length(which(out==1))>1) {
      omega2 = c(omega2, apply(modelos0.df[,which(out==1)], MARGIN = 2, 
                               function (x) ((x[4]*x[11] - x[2]*x[13]) -(x[4]+x[11])^2/4 )) )
    } else if (length(which(out==1))==1) {
      m1.aux = modelos0.df[,which(out==1)]   # aux vector with models coefs
      omega2 = c(omega2,(( m1.aux[4]*m1.aux[11] - m1.aux[2]*m1.aux[13]) 
                         -(m1.aux[4]+m1.aux[11])^2/4)  )
    }
    
  }
}
# Calculate period for all linear systems
period = data.frame(rep(años.ch,length(omega2)),  2*pi/sqrt(omega2))
colnames(period) = c("ENSO","T") 
# to obtain figure with EL Niño and La Niña distributions, run script for El Niño years 
# and then for La Niña years (do not overwrite period data.frame), and rbind new period
# with previous period 
# Usar el siguiente código para distribuciones ENSO
gg.samp = ggplot(data = period) + geom_histogram(aes(x = T, y=after_stat(density),
                                                     col = ENSO, fill = ENSO ),
                                                 binwidth = 1, position = "identity",
                                                 alpha = 0.4) +
  theme_bw() +
  theme(axis.text.x  = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 16, face="bold")) +
  theme(axis.text.y  = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size=16, face="bold")) +
  theme(legend.text = element_text(size = 15, face = "bold"), legend.direction = "horizontal",
        legend.position = "bottom",
        legend.title = element_text(size = 15, face = "bold")) +
  scale_fill_viridis_d(name = "", labels = c("La Niña","El Niño"),begin =0.1, end = 0.6,  direction = -1) +
  scale_color_viridis_d(name = "", labels = c("La Niña","El Niño"),begin = 0.1, end = 0.6, direction = -1) +
  xlab("Period [days]") + ylab("Frequency")
#
stats.0 = desc_statby(data = period, measure.var = "T", grps = "ENSO")
stats.samp = round(stats.0[,c("length","min","max", "mean","sd")]*10)/10
colnames(stats.samp)= c("#LM","Min","Max","Mean","Sd")
gg.table = ggtexttable(format(stats.samp, digits = 3), rows = NULL,
                       theme = ttheme(
                         colnames.style = colnames_style(color = "white", size = 14),
                         tbody.style = tbody_style(fill = viridis(2, begin=0.1, end = 0.6, direction = -1, alpha = 0.4), color = "black", size = 14),
                       ))

gg.samp + annotation_custom(ggplotGrob(gg.table), xmin = 67, ymin = 0.16)
# To save figure, uncomment the follwoing lines:
# ggsave("Period_enso.png", device = png(), height = 8, width = 10, 
#       dpi = 300, bg = "white")

# Plot structure of most frequent models:
modelos.ind.0 = which(freq.modelos*length(modelos0.df)>length(modelos0.df)*0.01)   # select models with fr>0.01
freq.sort = sort(freq.modelos[modelos.ind.0], decreasing = T, index.return = T)    # order models with decreasing frequency
modelos.ind = modelos.ind.0[freq.sort$ix]       # order rows of selected models with decreasing frequency
Modelos.df = cbind(seq(1,2*NPsi), modelos.df[,modelos.ind])
if(años.ch == "all") {
  colnames(Modelos.df) = c("Coef", paste0("C.M", seq(1,length(modelos.ind))))
} else if (años.ch == "niños") {
  colnames(Modelos.df) = c("Coef", paste0("EN.M", seq(1,length(modelos.ind))))
} else if (años.ch == "niñas") {
  colnames(Modelos.df) = c("Coef", paste0("LN.M", seq(1,length(modelos.ind))))
}

# Plot most frequently model's structure
Modelos.long.df = melt(data = Modelos.df, id.vars = "Coef") # change data.frame to long format
facet.labs = paste0(colnames(Modelos.df[-1]), 
                    rep(" fr=",length(modelos.ind)), round(freq.sort$x*100)/100 )
names(facet.labs) = paste0(colnames(Modelos.df[-1]))

gg.mask = ggplot() + geom_col(data = Modelos.long.df, aes(x = Coef, y = value, 
                                                               group = Coef, fill = ifelse(Coef<=NPsi,"fc1","fc2"))) +
  scale_fill_manual(values = c("fc1"="cornflowerblue", "fc2"="darkolivegreen3"), 
                    labels = c(TeX("$\\dot{x}$"),TeX("$\\dot{y}$")), name = "Equation" ) +
  scale_x_continuous(breaks = seq(1,18), labels = x.lab) +
  scale_y_continuous(breaks = c(0,1), labels = NULL ) + ylab(NULL) +
  theme(axis.text.x  = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.ticks = element_blank()) +
  theme(axis.text.y  = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size=16, face="bold")) +
  theme(legend.text = element_text(size = 15, face = "bold"), legend.direction = "horizontal",
        legend.position = "top", legend.spacing.x = unit(20,"pt"),
        legend.title = element_text(size = 15, face = "bold")) +
  facet_grid(variable~., switch = "y", labeller = labeller(variable = facet.labs)) + 
  theme(strip.text.y.left = element_text(angle = 0, size = 15),
        strip.text = element_text(margin = margin(t =10, r = 10, b = 10, l = 10)) )

# ggarrange(gg.mask.niño + rremove("xlab") + rremove("x.text") + rremove("x.ticks"),
#          gg.mask.niña, nrow = 2, align = "hv", common.legend = T, legend = "top", 
#          heights = c(0.99,1.1) ) + 
#  annotation_custom(grid.polygon(x = c(0,1), y = c(0.5,0.5), gp = gpar(col = "black", lty = "dashed", lwd = 2) ))

# ggsave("Modelos_masks.png", device = png(), height = 7, width =10, 
#       dpi = 300, bg = "white")
