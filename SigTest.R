### SigTest.R                               ###
### Author: Phil Woodward                   ###
### Pilot test code for Sams website        ###
### (last updated 2014 Jan 07)              ###
###############################################
### Reset R Environment ###
rm(list=ls())
ls()
###########################

#################################################################################
### SCENARIO SETTINGS (START)

DeltaLim = 0                       ### "sup" or "non-inf" delta limit for true delta. (must be >=0)
alpha = 0.05                       ### all 'tests' are 1-sided - set 1 value only

sigma = 10		      	     ### Residual st.dev. (between for PG or within for XO)
N     = 12
VarD = 2*sigma*sigma/N	           ### Enter as a function of sigma & N
#####################################################################################################
### Var.DeltaEst is the sampling variance of the conventional likelihood based estimator of delta ###
### Some typical values for common designs:									  ###
### PG with N arms per group: 2 res.var.between / N                                               ###
### XO with all N subjects receiving both treatments once only: 2 res.var.within / N              ###
### XO with all N subjects receiving one treatment once & the other twice: 1.5 res.var.within / N ###
### For more complex designs, Genstat can be used to obtain variance as constant times sigma^2    ###
#####################################################################################################
DF = 2*(N-1)		### Typically, PG= grps*(N-1), XO= prds*N-N-(trts-1)-(prds-1)

delta = seq(0, 25, by=1)  ### true value of delta (positive implies test is BETTER than comparator/control)
                          ### DeltaTrue will be on x-axis of graph, so should select many values to smooth the curve
### SCENARIO SETTINGS (END)
#################################################################################

###################################################
### Power calculations done using non-central t ###
Power = pt( qt(1-alpha, df= DF), df= DF, ncp= (delta - DeltaLim)/sqrt(VarD), lower.tail= FALSE )
###################################################

#############################################
###             PLOTS                     ###
###     ***** EDIT titles *****           ###
#############################################
### PLOT 1: Expected Decision Boundary    ###
### Normal approximation used here        ###
#############################################
### add later ###

########################################
###      Plot the power curve        ###
###      X-axis is delta.true        ###
### ***** EDIT titles & abline ***** ###
########################################
plot(delta, Power, type= "n", ylim= c(0,1), axes=F, xlab= "True Value of Delta", ylab= "Power")
title(main=paste("Power Curves (alpha:", alpha, ", H0 delta:", DeltaLim,")"), line=0, cex.main= 1)
axis(1)
axis(2, at= seq(0,1, by= 0.1))

lines(delta, Power, lwd= 2, lty= 1, col= "black")

grid()
abline(h=c(0.5,0.8), lty=3)
legend("bottomright", legend= paste( "N", N), lty= 1, bty= 'n', lwd= 2, col= "black")