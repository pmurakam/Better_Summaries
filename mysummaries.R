################################################################################
## Copyright (C) 2010 Peter Murakami <pmurakam@jhsph.edu>
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

#Function for better output from gee (package gee) function:
gee.mysummary <- function(fit, alpha=.05, dig=3, p.dig=4){
#fit is a fitted gee object. dig is number of digits to report.
    zq <- qnorm(1-alpha/2)
    estimate <- coef(fit)
    lower    <- coef(fit) - zq*summary(fit)$coef[,"Robust S.E."]
    upper    <- coef(fit) + zq*summary(fit)$coef[,"Robust S.E."]
    
    #robust.z  <- round(summary(fit)$coef[,"Robust z"], dig)
    robust.se <- round(summary(fit)$coef[,"Robust S.E."], dig)
    naive.se  <- round(summary(fit)$coef[,"Naive S.E."],dig)
    p <- 2*pnorm(abs(summary(fit)$coef[,"Robust z"]), lower.tail=FALSE)
    p <- signif(p, digits=p.dig)
    if(all(fit$family[1:2]==c("binomial","logit"))){
        Odds.Ratio <- round(exp(estimate), dig)
        OR.lower   <- round(exp(lower), dig)
        OR.upper   <- round(exp(upper), dig)
        estimate   <- round(estimate, dig)
        return(cbind(estimate,Odds.Ratio,OR.lower,OR.upper,
                     naive.se,robust.se,p))
    } else if(all(fit$family[1:2]==c("poisson","log"))){
        #assumes all counts are taken over same lengths of time 
        #or that an offset term has been specified, as usual.
        IRR       <- round(exp(estimate), dig) #incidence rate ratio
        IRR.lower <- round(exp(lower), dig)
        IRR.upper <- round(exp(upper), dig)
        estimate  <- round(estimate, dig)
        return(cbind(estimate,IRR,IRR.lower,IRR.upper,
                     naive.se,robust.se,p))
    } else{
        estimate <- round(estimate, dig)
        lower    <- round(lower, dig)
        upper    <- round(upper, dig)
        return(cbind(estimate,lower,upper,naive.se,robust.se,p))
    }
}

geese.mysummary <- function(fit, alpha=.05, dig=3, p.dig=4){
#fit is a fitted geese() object. dig is number of digits to report.
    zq <- qnorm(1-alpha/2)
    estimate <- fit$beta
    lower    <- fit$beta - zq*summary(fit)$mean[,"san.se"]
    upper    <- fit$beta + zq*summary(fit)$mean[,"san.se"]
    
    #robust.z  <- round(summary(fit)$coef[,"Robust z"], dig)
      #this robust.z is how gee() calculates it, but i exclude 
      #it here since i'm unsure what the Wald column means in 
      #the geese() output.
    robust.se <- round(summary(fit)$mean[,"san.se"], dig)
    p <- summary(fit)$mean[,"p"]
    p <- signif(p, digits=p.dig)
    if(all(fit$model[1:2]==c("logit","binomial"))){
        Odds.Ratio <- round(exp(estimate), dig)
        OR.lower   <- round(exp(lower), dig)
        OR.upper   <- round(exp(upper), dig)
        estimate   <- round(estimate, dig)
        return(cbind(estimate,robust.se,Odds.Ratio,OR.lower,OR.upper,p))
    } else if(all(fit$model[1:2]==c("log","poisson"))){
        #assumes all counts are taken over same lengths of time 
        #or that an offset term has been specified, as usual.
        IRR       <- round(exp(estimate), dig) #incidence rate ratio
        IRR.lower <- round(exp(lower), dig)
        IRR.upper <- round(exp(upper), dig)
        estimate  <- round(estimate, dig)
        return(cbind(estimate,robust.se,IRR,IRR.lower,IRR.upper,p))
    } else{
        estimate <- round(estimate, dig)
        lower    <- round(lower, dig)
        upper    <- round(upper, dig)
        return(cbind(estimate,robust.se,lower,upper,p))
    }
}

#Function for better glm output:
#see p.31 in Extending the Linear Model with R by Faraway for 
#difference between confidence interval types.
glm.mysummary <- function(fit, alpha=.05, dig=3, ci="profile.lik", p.dig=4){
#fit is a fitted gee object. dig is number of digits to report.
    zq <- qnorm(1-alpha/2)
    estimate <- coef(fit)
    #lower    <- coef(fit) - zq*summary(fit)$coef[,"Std. Error"]
    #upper    <- coef(fit) + zq*summary(fit)$coef[,"Std. Error"]
    if(ci=="profile.lik"){
        require(MASS)
        cin <- confint(fit)
        lower <- cin[,1]
        upper <- cin[,2]
    } else if(ci=="normal.approx"){ 
        lower <- coef(fit) - zq*summary(fit)$coef[,"Std. Error"]
        upper <- coef(fit) + zq*summary(fit)$coef[,"Std. Error"] 
    }
    
    se <- round(summary(fit)$coef[,"Std. Error"], dig)
    st <- which(colnames(summary(fit)$coef)%in%c("z value","t value"))
    st <- colnames(summary(fit)$coef)[st]
    stats <- round(summary(fit)$coef[,st], dig)
    pcol <- grep("Pr",colnames(summary(fit)$coef))
    p  <- summary(fit)$coef[,pcol]
    p  <- signif(p, digits=p.dig)
    if(all(fit$family[1:2]==c("binomial","logit"))){
        Odds.Ratio <- round(exp(estimate), dig)
        OR.lower   <- round(exp(lower), dig)
        OR.upper   <- round(exp(upper), dig)
        estimate   <- round(estimate, dig)
        output <- cbind(estimate,Odds.Ratio,OR.lower,OR.upper,stats,p)
        colnames(output)[5] <- st
        return(output)
    } else if(all(fit$family[1:2]==c("poisson","log"))){
        #assumes all counts are taken over same lengths of time 
        #or that an offset term has been specified, as usual.
        IRR       <- round(exp(estimate), dig) #incidence rate ratio
        IRR.lower <- round(exp(lower), dig)
        IRR.upper <- round(exp(upper), dig)
        estimate  <- round(estimate, dig)
        output <- cbind(estimate,IRR,IRR.lower,IRR.upper,stats,p)
        colnames(output)[5] <- st
        return(output)
    } else{
        estimate <- round(estimate, dig)
        lower    <- round(lower, dig)
        upper    <- round(upper, dig)
        output <- cbind(estimate,lower,upper,stats,p)
        colnames(output)[4] <- st
        return(output)
    }
}

lm.mysummary <- function(fit, alpha=.05, dig=3, p.dig=4){
#fit is a fitted gee object. dig is number of digits to report.
    zq <- qnorm(1-alpha/2)
    estimate <- round(coef(fit), dig)
    lower <- round(coef(fit)- zq*summary(fit)$coef[,"Std. Error"], dig)
    upper <- round(coef(fit)+ zq*summary(fit)$coef[,"Std. Error"], dig)
    t  <- round(summary(fit)$coef[,"t value"], dig)
    se <- round(summary(fit)$coef[,"Std. Error"], dig)
    p  <- signif(summary(fit)$coef[,"Pr(>|t|)"], digits=p.dig)
    
    return(cbind(estimate,se,lower,upper,t,p))
}
