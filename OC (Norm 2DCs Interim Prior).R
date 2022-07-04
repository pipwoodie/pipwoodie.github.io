### OC (Norm 2DCs Interim Prior).R                   ###
### Author: Phil Woodward                            ###
### (last updated 2013 Dec 08)                       ###
### Operating Characteristics for Normal case        ###
### with 2 Decision Criteria,                        ###
### Bayesian Interim Analysis & informative prior    ###
### x-axis= true delta                               ###
### Modified for ??? Project by ???, 20?? ??? ??     ###
###########################################################################
### Ref: Walley, R.J., Smith, C.L., Gale, J.D. and Woodward, P. (2015). ###
### Advantages of a wholly Bayesian approach to assessing efficacy in   ###
### early drug development: a case study.                               ###
### Pharmaceutical Statistics, Vol 14, Iss 3, pp205-215.                ###
###########################################################################
### Reset R Environment ###
rm(list=ls())
ls()

#################################################################################
### SCENARIO SETTINGS (START) - all settings are scalar unless stated otherwise
### For criterion, at least 50% sure True Delta greater than 2: DeltaLim1= 2 & Prob1= 0.50 ("Mean level setting")
DeltaLim1 = 4                     ### First delta limit for true delta. (must be >=0)
Prob1 = 0.50                      ### Probability level associated with DeltaLim1.
### For criterion, at least 75% sure True Delta greater than 1: DeltaLim2= 1 & Prob2= 0.75 ("Precision level setting")
DeltaLim2 = 2                    ### Second delta limit for true delta. (must be >=0) - set 1 value only
Prob2 = 0.75                     ### Probability level associated with DeltaLim2 - set 1 value only
### A Bayesian interpretation of success criteria: Pr[ True Delta > DeltaLim# ] > Prob# ###
###################
### Study should be designed so that if ('mean level') Criterion 1 is met, then ('precision level') Criterion 2 will also.
### Hence, Criterion 1 will dominate the OC and this should be the focus when assessing the design
###################

### Interim Predictive Probability Thresholds
PredProbLim.futile  = 0.2   ### = 0 if only considering success at interim (very unlikely)
PredProbLim.success = 0.8   ### = 1 if only considering futility at interim
 
sigma = 5
Ntot.Ctrl = 25         ### Number subjects in study:
Ntot.Actv = 50         ### Ntot = total at end of study
Nint.Ctrl = 12         ### Nint = total at interim,  Npt2 = total recruited after the interim
Nint.Actv = 25         ### .Ctrl = Control group,  .Actv = Active group
### DO NOT EDIT FOLLOWING FORMULAE ###
Npt2.Ctrl = Ntot.Ctrl - Nint.Ctrl  ###
Npt2.Actv = Ntot.Actv - Nint.Actv  ###
######################################
 
### Prior Distributions for Control Mean Response ###
### Vague prior can be specified by setting PriorA.CtrlSD = 1000*sigma ###
PriorA.CtrlMean = 0        ### PriorA = Analysis prior that will be used in formal analysis
PriorA.CtrlSD   = 1        ### PriorD = Design prior that will only be used in design assessment calculations
PriorD.CtrlMean = PriorA.CtrlMean  ### .CtrlMean = Prior Mean for Control Mean Response
PriorD.CtrlSD = PriorA.CtrlSD      ### .CtrlSD = Prior St.Dev. for Control Mean Response
################################################################################################
### Only have different Design Prior if want to assess impact of mis-specifed Analysis Prior ###
### DO NOT EDIT FOLLOWING FORMULAE ##############################
PriorA.CtrlEffN = (sigma / PriorA.CtrlSD)^2                   ###
PriorD.CtrlEffN = (sigma / PriorD.CtrlSD)^2                   ###
#################################################################
 
v.delta = seq(0, 8, by=0.25)  ### Vector giving range of true values of delta (+ve implies Active is BETTER than Control)
                              ### delta will be on x-axis of graph, so should select many values to smooth the curve
                              ### However, this will increase the number of simulations being run, so recommend 20 to 40.
 
N.sims <- 10000              ### No. simulations used to estimate OC Curve (1000 at first, 10000 once design fixed)
 
### SCENARIO SETTINGS (END)
### Any Warnings? ###
if ( min(v.delta) < 0 ) { print( "WARNING: a negative DeltaTrue.val implies test treatment is worse than comparator/control" ) }
if ( Npt2.Ctrl < 1 ) { print( "ERROR: Npt2.Ctrl must be >0; CALCULATIONS WILL BE INCORRECT!" ) }
if ( Npt2.Actv < 1 ) { print( "ERROR: Npt2.Actv must be >0; CALCULATIONS WILL BE INCORRECT!" ) }
#################################################################################
 
delta.len = length(v.delta)
m.delta   = matrix( rep(v.delta, N.sims), nrow= N.sims, byrow=T)
endPostSD = sigma * sqrt( (PriorA.CtrlEffN + Ntot.Ctrl + Ntot.Actv) / 
                            (Ntot.Actv * (PriorA.CtrlEffN + Ntot.Ctrl)) )
intPostSD = sigma * sqrt( (PriorA.CtrlEffN + Nint.Ctrl + Nint.Actv) / 
                            (Nint.Actv * (PriorA.CtrlEffN + Nint.Ctrl)) )
 
### Simulate sample statistics at the interim.  First need Control Mean from Design Prior. ###
m.CtrlMean <- matrix(rnorm(N.sims*delta.len, mean= PriorD.CtrlMean, sd= PriorD.CtrlSD), nrow= N.sims, byrow= T)
m.CtrlXbar <- matrix(rnorm(N.sims*delta.len, mean= 0,       sd= sigma/sqrt(Nint.Ctrl)), nrow= N.sims, byrow= T) + m.CtrlMean
m.ActvXbar <- matrix(rnorm(N.sims*delta.len, mean= v.delta, sd= sigma/sqrt(Nint.Actv)), nrow= N.sims, byrow= T) + m.CtrlMean

### Predictive distribution at interim for linear combination of post interim sample means
m.intPredMean = (Npt2.Actv / Ntot.Actv) * m.ActvXbar -
                (Npt2.Ctrl /(PriorA.CtrlEffN + Ntot.Ctrl)) *
                ((PriorA.CtrlEffN * PriorA.CtrlMean + Nint.Ctrl * m.CtrlXbar) /
                 (PriorA.CtrlEffN + Nint.Ctrl) )

intPredVar = ( (Npt2.Actv / Ntot.Actv)^2 * ((1/Nint.Actv) + (1/Npt2.Actv)) +
                 (Npt2.Ctrl /(PriorA.CtrlEffN + Ntot.Ctrl))^2 *
                 ((1/(PriorA.CtrlEffN + Nint.Ctrl)) + (1/Npt2.Ctrl)) ) * sigma^2

intPredSD = sqrt(intPredVar)

### Evaluate RHS of inequality, and then Predictive Probability of end of study success
### Criterion 1 ###
m.RHS1 = DeltaLim1 + qnorm(Prob1) * endPostSD + 
         (PriorA.CtrlEffN * PriorA.CtrlMean + Nint.Ctrl * m.CtrlXbar) / (PriorA.CtrlEffN + Ntot.Ctrl) -
         (Nint.Actv / Ntot.Actv) * m.ActvXbar

m.PredProbSuccess1 = pnorm( m.RHS1, mean= m.intPredMean, sd= intPredSD, lower.tail=F)
m.intSuccess1 = (m.PredProbSuccess1 > PredProbLim.success)
m.intFutile1  = (m.PredProbSuccess1 < PredProbLim.futile)
m.intCont1    = !(m.intSuccess1 | m.intFutile1)

v.intProbSuccess1 = apply(m.intSuccess1, 2, mean)
v.intProbFutile1  = apply(m.intFutile1, 2, mean)
v.intProbCont1    = 1 - v.intProbSuccess1 - v.intProbFutile1
v.intProbStop1    = 1 - v.intProbCont1
######################################
### Criterion 2 ###
m.RHS2 = DeltaLim2 + qnorm(Prob2) * endPostSD + 
         (PriorA.CtrlEffN * PriorA.CtrlMean + Nint.Ctrl * m.CtrlXbar) / (PriorA.CtrlEffN + Ntot.Ctrl) -
         (Nint.Actv / Ntot.Actv) * m.ActvXbar

m.PredProbSuccess2 = pnorm( m.RHS2, mean= m.intPredMean, sd= intPredSD, lower.tail=F)
m.intSuccess2 = (m.PredProbSuccess2 > PredProbLim.success)
m.intFutile2  = (m.PredProbSuccess2 < PredProbLim.futile)
m.intCont2    = !(m.intSuccess2 | m.intFutile2)

v.intProbSuccess2 = apply(m.intSuccess2, 2, mean)
v.intProbFutile2  = apply(m.intFutile2, 2, mean)
v.intProbCont2    = 1 - v.intProbSuccess2 - v.intProbFutile2
v.intProbStop2    = 1 - v.intProbCont2
######################################
### Both Criteria ###
### Approximate by noting one criterion will dominate for known sigma case
m.intSuccessB = m.intSuccess1 & m.intSuccess2  ## both must succeed
m.intFutileB  = m.intFutile1  | m.intFutile2   ## either can be futile
m.intContB    = !(m.intSuccessB | m.intFutileB)

v.intProbSuccessB = apply(m.intSuccessB, 2, mean)
v.intProbFutileB  = apply(m.intFutileB, 2, mean)
v.intProbContB    = 1 - v.intProbSuccessB - v.intProbFutileB
v.intProbStopB    = 1 - v.intProbContB
######################################

### Evaluate probability of success at interim conditional on delta; needed for conventional OC curve
### Averaging over prior for CtrlMean using simulated values generated earlier

m.intCondMean = (Npt2.Actv / Ntot.Actv) * (m.CtrlMean + m.delta) -
                (Npt2.Ctrl /(PriorA.CtrlEffN + Ntot.Ctrl)) * m.CtrlMean

intCondVar = ( (Npt2.Actv / Ntot.Actv^2) + (Npt2.Ctrl /(PriorA.CtrlEffN + Ntot.Ctrl)^2)  ) * sigma^2

intCondSD = sqrt(intCondVar)

m.ConProbSuccess1 = pnorm( m.RHS1, mean= m.intCondMean, sd= intCondSD, lower.tail=F)
m.ConProbSuccess2 = pnorm( m.RHS2, mean= m.intCondMean, sd= intCondSD, lower.tail=F)
### Probability both criteria met approximated by lowest of either as known sigma assumption causes one to dominate
m.ConProbSuccessB = min(m.ConProbSuccess1, m.ConProbSuccess2)

### For comparison, calculate the Power assuming no interim analysis ###
v.PriorPredDeltaMean.Mean = PriorD.CtrlMean + v.delta - 
                         (PriorA.CtrlEffN * PriorA.CtrlMean + Ntot.Ctrl * PriorD.CtrlMean) / (PriorA.CtrlEffN + Ntot.Ctrl)
v.PriorPredDeltaMean.SD   = sqrt( (1/Ntot.Actv) + (Ntot.Ctrl / (PriorA.CtrlEffN + Ntot.Ctrl)^2) +
                                  (PriorA.CtrlEffN / (PriorA.CtrlEffN + Ntot.Ctrl))^2 / PriorD.CtrlEffN ) * sigma

v.Power1.NoInt = pnorm(DeltaLim1 + qnorm(Prob1) * endPostSD, mean= v.PriorPredDeltaMean.Mean, 
                      sd= v.PriorPredDeltaMean.SD, lower.tail= F)
v.Power2.NoInt = pnorm(DeltaLim2 + qnorm(Prob2) * endPostSD, mean= v.PriorPredDeltaMean.Mean, 
                      sd= v.PriorPredDeltaMean.SD, lower.tail= F)

### Following is same as conditional probability of success at study start assuming no interim analysis
### Only use is to check simulation code is correct
##v.Power1.NoInt.sim = apply(m.ConProbSuccess1, 2, mean)
##v.Power2.NoInt.sim = apply(m.ConProbSuccess2, 2, mean)
### Visual/crude check simulation is good approximation to exact result when no interim ###
##plot(v.Power1.NoInt.sim, v.Power1.NoInt)
##100*(v.Power1.NoInt.sim - v.Power1.NoInt)
##plot(v.Power2.NoInt.sim, v.Power2.NoInt)
##100*(v.Power2.NoInt.sim - v.Power2.NoInt)
###########################################################################################

### ProbSuccess, conditional on known delta & interim decision
### ProbNotFutile x ProbSuccess | NotFutile
m.ConProbSuccess1.WithInt = m.ConProbSuccess1 * m.intCont1 + m.intSuccess1
v.ConProbSuccess1.WithInt   = apply(m.ConProbSuccess1.WithInt, 2, mean)
m.ConProbSuccess2.WithInt = m.ConProbSuccess2 * m.intCont2 + m.intSuccess2
v.ConProbSuccess2.WithInt   = apply(m.ConProbSuccess2.WithInt, 2, mean)
m.ConProbSuccessB.WithInt = m.ConProbSuccessB * m.intContB + m.intSuccessB
v.ConProbSuccessB.WithInt   = apply(m.ConProbSuccessB.WithInt, 2, mean)
##########################################################################

########################################
###             PLOTS                ###
###      X-axis is delta.true        ###
###     ***** EDIT titles *****      ###
########################################
###    PLOT 1: OC at Interim         ###
########################################
par(mfrow=c(2,2))
### Graph (1,1) ###
plot(c(0,1), c(0,1), type= "n", axes=F, xlab="", ylab="")
title(main= paste("OC at INTERIM"), line=0, cex.main= 1)
text(x=0,y=0.95,pos=4, labels= paste("Pred Probs defining Futility/Success:", 
                                    PredProbLim.futile, "/", PredProbLim.success))
text(x=0,y=0.85,pos=4, labels= paste("Criterion 1: Prob > ", Prob1, " that delta > ", DeltaLim1, " (solid line)", sep=""))
text(x=0,y=0.75,pos=4, labels= paste("Criterion 2: Prob > ", Prob2, " that delta > ", DeltaLim2, " (dashed line)", sep=""))
text(x=0,y=0.65,pos=4, labels= "Both criteria: broken red line (may overlap Criterion 1 or 2)")
text(x=0,y=0.5,pos=4, labels= paste("Interim N.  Control: ", Nint.Ctrl, " of ", Ntot.Ctrl, 
                                    ",  Active: ", Nint.Actv, " of ", Ntot.Actv, sep=""))
text(x=0,y=0.35,pos=4, labels= paste( "Prior for Control Mean:  N( mean= ", PriorA.CtrlMean, 
                                                                  ", sd= " , PriorA.CtrlSD, " )", sep=""))
text(x=0,y=0.25,pos=4, labels= paste( "Prior Effective N (Controls) = ", format(PriorA.CtrlEffN, digits=1), sep=""))
text(x=0,y=0.1,pos=4, labels= paste( "sigma = ", sigma, ", Effect's PostSD@int = ", format(intPostSD, digits=2), sep=""))

### Graph (2,1) ###
plot(v.delta, v.intProbSuccess1, type= "n", ylim= c(0,1), axes=F,     
     xlab= "True Value of Delta", ylab= "Probability Success")
title(main= "SUCCESS", line=0, cex.main= 1)
axis(1)
axis(2, at= seq(0,1, by= 0.1))
lines(v.delta, v.intProbSuccess1, lty= 1, lwd= 2, col= 1)
lines(v.delta, v.intProbSuccess2, lty= 2, lwd= 2, col= 1)
lines(v.delta, v.intProbSuccessB, lty= 3, lwd= 2, col= "red")
grid()
abline(h=c(0.5,0.8), v=c(DeltaLim1, DeltaLim2), lty=3)

### Graph (1,2) ###
plot(v.delta, v.intProbCont1, type= "n", ylim= c(0,1), axes=F,     
     xlab= "True Value of Delta", ylab= "Probability Continue")
title(main= "CONTINUE", line=0, cex.main= 1)
axis(1)
axis(2, at= seq(0,1, by= 0.1))
lines(v.delta, v.intProbCont1, lty= 1, lwd= 2, col= 1)
lines(v.delta, v.intProbCont2, lty= 2, lwd= 2, col= 1)
lines(v.delta, v.intProbContB, lty= 2, lwd= 2, col = "red")
grid()
abline(h=c(0.5,0.8), v=c(DeltaLim1, DeltaLim2), lty=3)

### Graph (2,2) ###
plot(v.delta, v.intProbFutile1, type= "n", ylim= c(0,1), axes=F,     
     xlab= "True Value of Delta", ylab= "Probability Futile")
title(main= "FUTILITY", line=0, cex.main= 1)
axis(1)
axis(2, at= seq(0,1, by= 0.1))
lines(v.delta, v.intProbFutile1, lty= 1, lwd= 2, col= 1)
lines(v.delta, v.intProbFutile2, lty= 2, lwd= 2, col= 1)
lines(v.delta, v.intProbFutileB, lty= 2, lwd= 2, col= "red")
grid()
abline(h=c(0.5,0.8), v=c(DeltaLim1, DeltaLim2), lty=3)

########################################
###    PLOT 2: OC at Study End       ###
### 2 separate plots                 ###
########################################
par(mfrow=c(2,1))
### Graph (1,1) ###
plot(v.delta, v.ConProbSuccess1.WithInt, type= "n", ylim= c(0,1), axes=F,     
     xlab= "True Value of Delta", ylab= "Probability Go Decision")
title(main= paste("OC with Interim. Futility/Success Pred Prob Limits: ", PredProbLim.futile, "/", PredProbLim.success, sep=""))
axis(1)
axis(2, at= seq(0,1, by= 0.1))
lines(v.delta, v.ConProbSuccess1.WithInt, lty= 1, lwd= 2, col= 1)
lines(v.delta, v.ConProbSuccess2.WithInt, lty= 2, lwd= 2, col= 1)
grid()
abline(h=c(0.5,0.8), v=c(DeltaLim1, DeltaLim2), lty=3)

### Graph (2,1) ###
plot(v.delta, v.Power1.NoInt, type= "n", ylim= c(0,1), axes=F,     
     xlab= "True Value of Delta", ylab= "Probability Go Decision")
title(main= paste("OC without Interim. Criterion1(P>", Prob1, ", D>", DeltaLim1, ") Criterion2(P>", Prob2, ", D>", DeltaLim2, ")", sep=""))
axis(1)
axis(2, at= seq(0,1, by= 0.1))
lines(v.delta, v.Power1.NoInt, lty= 1, lwd= 2, col= 1)
lines(v.delta, v.Power2.NoInt, lty= 2, lwd= 2, col= 1)
grid()
abline(h=c(0.5,0.8), v=c(DeltaLim1, DeltaLim2), lty=3)

legend("bottomright", legend= c( paste( "sigma (PostSD): ", sigma,
                                 " (", format(endPostSD , digits=2), ")", sep="" ),
                                 paste("C1:",Prob1,">",DeltaLim1), paste("C2",Prob2,">",DeltaLim2)), lty = c(0,1,2), lwd= 2, bty= 'n')

########################################
###    PLOT 2: OC at Study End       ###
### both on same plot                ###
########################################
plot(v.delta, v.ConProbSuccess1.WithInt, type= "n", ylim= c(0,1), axes=F,     
     xlab= "True Value of Delta", ylab= "Probability Go Decision")
title(main= paste("OC with/without Interim. Futility/Success Pred Prob Limits: ", PredProbLim.futile, "/", PredProbLim.success, sep=""))
axis(1)
axis(2, at= seq(0,1, by= 0.1))
lines(v.delta, v.ConProbSuccess1.WithInt, lty= 1, lwd= 2, col= "blue")
lines(v.delta, v.ConProbSuccess2.WithInt, lty= 2, lwd= 2, col= "blue")
lines(v.delta, v.Power1.NoInt, lty= 1, lwd= 2, col= "red")
lines(v.delta, v.Power2.NoInt, lty= 2, lwd= 2, col= "red")
grid()
abline(h=c(0.5,0.8), v=c(DeltaLim1, DeltaLim2), lty=3)

grid()
abline(h=c(0.5,0.8), v=c(DeltaLim1, DeltaLim2), lty=3)

legend("bottomright", legend= c( paste( "sigma (PostSD): ", sigma,
                                 " (", format(endPostSD , digits=2), ")", sep="" ),
                                 paste("C1:",Prob1,">",DeltaLim1), paste("C2",Prob2,">",DeltaLim2), "with Interim", "without Interim"), 
                              lty = c(0,1,2,1,1), lwd= 2, bty= 'n', col= c("black","black","black","blue","red") )

#############################################
### PLOT 3: Expected Decision Boundaries  ###
###         at study end                  ###
#############################################
C1.mean = DeltaLim1 + qnorm(Prob1)*endPostSD
C2.mean   = DeltaLim2 + qnorm(Prob2)*endPostSD
lb.Xaxis = "Curves represent Probability Distribution of Effect (Posterior to Study)"

plot(y= c(1,5), x= c(min(v.delta), max(v.delta)), type='n', axes= F,
     ylab=paste("sigma:", sigma, "    N(Ctrl/Actv):", Ntot.Ctrl, "/", Ntot.Actv,
                " Prior Eff.N(Ctrl):", format(PriorA.CtrlEffN, digits=1),
                "    PostSD:", format(endPostSD, digits=2)), 
     xlab=lb.Xaxis, main= "Decision Criteria: Minimum Evidence Required to meet Decision Criteria")
axis(2, at= c(1,2,4,5), labels= c('', paste("C1:",Prob1,">",DeltaLim1), paste("C2",Prob2,">",DeltaLim2), ''))
axis(1)
abline(v= c(DeltaLim1, DeltaLim2), lty= 3)
lines(y= c(0,2), x= c(C1.mean, C1.mean), lty= 2, col= "blue")
lines(y= c(0,4), x= c(C2.mean, C2.mean), lty= 2, col= "dark green")

Vscale = 0.7/dnorm(0)

range1 = C1.mean + seq(-3, 3, by= 0.1)*endPostSD
pdf1 = 2 + dnorm(seq(-3, 3, by= 0.1))*Vscale
lines( range1, pdf1, col= "blue")
lines( c(min(range1),max(range1)), c(2,2), col= "blue")
lines( c(C1.mean - qnorm(Prob1)*endPostSD, C1.mean - qnorm(Prob1)*endPostSD), 
       c(2, 2 + dnorm(qnorm(Prob1))*Vscale), col= "blue")
points(C1.mean, 2, pch=19, col= "blue", cex= 3)

range2 = C2.mean + seq(-3, 3, by= 0.1)*endPostSD
pdf2 = 4 + dnorm(seq(-3, 3, by= 0.1))*Vscale
lines( range2, pdf2, col= "dark green")
lines( c(min(range2),max(range2)), c(4,4), col= "dark green")
lines( c(C2.mean - qnorm(Prob2)*endPostSD, C2.mean - qnorm(Prob2)*endPostSD), 
       c(4, 4 + dnorm(qnorm(Prob2))*Vscale), col= "dark green")
points(C2.mean, 4, pch=19, col= "dark green", cex= 3)

