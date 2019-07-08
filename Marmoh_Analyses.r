#Set working directory

setwd("~path-to-folder")

#load packages

library(survival)
library(survminer)
library(interval)
library(lme4)
library(brglm2)
library(glmmTMB)
library(coxme)
library(sandwich)
library(emmeans)
library(multcomp)
library(car)				
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)

#Set right contrasts for type III ss. 

options(contrasts = c("contr.sum", "contr.poly"))

########################################################################################
############################# Analysis 1 - time to event ###############################
########################################################################################

dat = read.csv("Marmoh_Analysis 1_T2E.csv", header = T)

#expand data set to have record for each seed in each dish at each sampling event. 
dat.expand = dat %>% slice(rep(1:n(),each=20))


success = rep(0,nrow(dat.expand))	#Blank vector in which to put binary success/failure
for (i in 1:150){					#initiate for() loop for each plate for each sample time (5 times 30)
nSuccess = dat[i,]$tGrm				#number of successes is equal to the total number of germinants for that dish for that week
nUnvia = dat[i,]$tUnv				#number of unviable seeds
ndChk = dat[i,]$dChk				#number of seeds with no event
if(i == 1){							#Begin if statement to modify indexing so that each row in original data frame puts its 20 values in the correct place. 
success[i:20] = c(rep(1,nSuccess),rep(0,ndChk),rep(2,nUnvia))	#Give as many 1-values as there are success, everything else is 0. 
} else { 
i2 = ((i-1)*20)+1					#modified index. At i=2, this starts at 21, at i=3 this starts at 41
success[i2:(i2+19)] = c(rep(1,nSuccess),rep(0,ndChk),rep(2,nUnvia)) #Same as with the first part of if statement, but places vector from i2 to 19 places after that.
}
}

dat.expand$Germinated = success	#Germinated = success
	
dat.expand2 = dat.expand

seedNum = rep(1:20,150)				#add seed number in dish for each dish observation time

dat.expand2$seedNum = paste0(dat.expand$dLabel,seedNum)	#create unique seed ID by combining dish name with number in dish

dat.expand2 = arrange(dat.expand2,seedNum)					#arrange observations by seed
dat.expand2.0 = dat.expand2[dat.expand2$Germinated == 0,]		#split data into obersavtions that had not germinated
dat.expand2.1 = dat.expand2[dat.expand2$Germinated != 0,]		#split data into observations that have germinated
dat.expand2.1 = arrange(dat.expand2.1,seedNum,Start)			#arrange germinated ones by seed ID again,and also Start time
dat.expand2.1x = dat.expand2.1[!duplicated(dat.expand2.1$seedNum),]	#Remove duplicated values, the first one (point of observed germination) is the only one kept

dat.expand2 = rbind(dat.expand2.0,dat.expand2.1x)	#combine the split datasets and overwrite old one
dat.expand2 = arrange(dat.expand2,seedNum,Start)	#arrange by seed number again. 

dat.expand2.germ = dat.expand2[dat.expand2$Germinated == 1,]								#Now keep only the observation at which it germinated or else...
dat.expand2.cens = dat.expand2[dat.expand2$Germinated == 0 & dat.expand2$End == 5,] #...keep the observation at which it became right-censored
dat.expand2.2 = rbind(dat.expand2.germ,dat.expand2.cens)	
dat.expand2.2.1 = rbind(dat.expand2.germ,dat.expand2.cens)		#Now combine them, and you have each observation, the time it germinated, or whether it never germinated
dat.expand2.2$End[dat.expand2.2$End == 5 & dat.expand2.2$Germinated == 0] <- Inf
dat.expand2.2$Start[dat.expand2.2$Start == 4 & dat.expand2.2$Germinated ==0] <-5
levels(dat.expand2.2$StratifLength) <- c("Nonstratified", "Cold Stratified (12 Weeks)")

#Can look at nonparametric curves
so.r = with(dat.expand2.2.1,Surv(End,event=Germinated,type="right")) #Right censored to use coxph to test proportional hazards assumption
so = with(dat.expand2.2,Surv(time = Start, time2 = End, type ="interval2"))
sf = survfit(so ~ IncubationTemp, data = dat.expand2.2)
ggsurvplot(sf,fun="event",xlim=c(0,5),xlab = "Weeks", ylab = "Cumulative Germination",legend.title="Treatment")

#Kaplan-Meier Curves with group assignments from log-rank tests (seen later)
surv_fit(so ~ StratifLength + IncubationTemp, dat.expand2.2) %>%
	surv_summary(dat.expand2.2) %>%
	add_row(time = 0,
		n.risk = tapply(.$n.risk, .$strata, max),
		n.event = 0, n.censor = 0,
		surv = 1, std.err = 0,
		upper = 1, lower = 1,
		strata = levels(.$strata),
		StratifLength = rep(levels(.$StratifLength),each = nlevels(.$IncubationTemp)),
		IncubationTemp = rep(levels(.$IncubationTemp), nlevels(.$StratifLength)), .before=T) %>%
	arrange(time,StratifLength,IncubationTemp) %>%
	mutate(Group = as.character(c(rep(NA,23),"a","d","a","b","d","c")))%>%
ggplot(aes(time, 1-surv, linetype = IncubationTemp,label=Group)) + theme_classic() + geom_step(color = "Black", size = 1.2) + 
facet_wrap(aes(group = StratifLength))+
scale_y_continuous("Cumulative Germination", labels = scales::percent)+
xlab("Weeks")+
guides(linetype=guide_legend(title=expression("Incubation Temperature " (degree~C))))+
scale_linetype_manual(values=c("solid","dotted","twodash"))+
geom_text(position=position_nudge(x=0.2))+
theme(panel.grid = element_blank(), panel.border = element_blank(), text = element_text(size = 18), legend.text=element_text(size=16),legend.title=element_text(size=16),
	axis.line=element_line(size=1),axis.text.y=element_text(size=16),axis.text.x=element_text(size=16),axis.title=element_text(size=18),
	plot.title = element_text(size=16), legend.position = "top")
ggsave("Marmoh_Exp1_KM_STEP.tiff",units="in",width=12,height=6,dpi=300)

#Test assumptions for parametric T2E models
sf2 = survfit(so ~ StratifLength + IncubationTemp, data = dat.expand2.2,type = "fleming-harrington")
time.sf = sf2$time
surv.sf=sf2$surv
trt = c(rep("StratifLength=0 Weeks, IncubationTemp=15/6",3),rep("StratifLength=0 Weeks, IncubationTemp=25/15",3),rep("StratifLength=0 Weeks, IncubationTemp=35/20",5),
rep("StratifLength=12 Weeks, IncubationTemp=15/6",5),rep("StratifLength=12 Weeks, IncubationTemp=25/15",4),rep("StratifLength=12 Weeks, IncubationTemp=35/20",3))
plotdat = data.frame(surv.sf = surv.sf, time.sf = time.sf, trt = trt)
ggplot(data = plotdat, aes(x=log(time.sf),y=log(-log(surv.sf))),group = trt) + theme_classic()+
geom_line(aes(color = trt))		#Weibull distribution does not fit, lines should be straight and not cross (i.e. proportional odds)

ggplot(data = plotdat, aes(x=log(time.sf),y=log(exp(-log(surv.sf))-1)),group = trt) + theme_classic()+
geom_line(aes(color = trt))		#similarly poor fit for log-logistic distribution


#First lets try semi-parametric Cox model
expt1.cox = coxph(so.r~StratifLength*IncubationTemp,data = dat.expand2.2) #build cox proportional hazards model
cox.zph(expt1.cox)																#Assumptions do not hold, many p-values below 0.05
ggcoxdiagnostics(expt1.cox,type = "deviance",linear.predictions=F)			#Diagnostic plot shows strong patterns in variance, not symmetric around 0
ggcoxdiagnostics(expt1.cox,type = "dfbeta",linear.predictions=F)			#Some observations are far more influential than others. 

#Parametric AFT models fit poorly, proportional odds and proportional hazards assumptions are not met, models shown below

ND = with(dat.expand2.2,expand.grid(StratifLength=levels(StratifLength),IncubationTemp=levels(IncubationTemp),tUnv=unique(tUnv)))

expt1.cox.pred = predict(expt1.cox, type = "risk") #Loglikelihood = -2222.219...this is very high compared to other models .

expt1.weib = survreg(so~StratifLength*IncubationTemp + frailty(dLabel),data=dat.expand2.2,dist="weibull") #over all is sig (Chi = 27.77, df = 5, p = 4e-5), LL = -724.9254
expt1.weib.aov = Anova(expt1.weib,test = "LR",type = 3) #Everything Sig, used type 3 b/c very unbalanced after censoring and want to examine interaction

expt1.llog = survreg(so~StratifLength * IncubationTemp + frailty(dLabel),data=dat.expand2.2,dist="loglogistic") #over all is sig (Chi = 27.77, df = 5, p = 4e-5), LL = -630.9722
expt1.llog.aov = Anova(expt1.llog,type = 3,test = "LR") #Everything Sig, used type 3 b/c very unbalanced after censoring and want to examine interaction

expt1.llog.pred = survreg(so~StratifLength*IncubationTemp,data=dat.expand2.2,dist="loglogistic")
#Data for plot

#Generate predictions and make plots for loglinear AFT Models
pt = predict(expt1.llog.pred, newdata = ND, type = "quantile", se=T, p = seq(0.0,0.999,0.001))
pt1 = data.frame(StratifLength = "0 Weeks", IncubationTemp = "15/6", time = pt$fit[1,],se=pt$se[1,],pct = 0:999/1000)
pt2 = data.frame(StratifLength = "12 Weeks", IncubationTemp = "15/6", time = pt$fit[2,],se=pt$se[2,],pct = 0:999/1000)
pt3 = data.frame(StratifLength = "0 Weeks", IncubationTemp = "25/15", time = pt$fit[3,],se=pt$se[3,],pct = 0:999/1000)
pt4 = data.frame(StratifLength = "12 Weeks", IncubationTemp = "25/15", time = pt$fit[4,],se=pt$se[4,],pct = 0:999/1000)
pt5 = data.frame(StratifLength = "0 Weeks", IncubationTemp = "35/20", time = pt$fit[5,],se=pt$se[5,],pct = 0:999/1000)
pt6 = data.frame(StratifLength = "12 Weeks", IncubationTemp = "35/20", time = pt$fit[6,],se=pt$se[6,],pct = 0:999/1000)
pt.full = rbind(pt1,pt2,pt3,pt4,pt5,pt6)
pt.full$up = pt.full$pct + pt.full$se
pt.full$lo = pt.full$pct - pt.full$se
pt.full$lo[pt.full$lo < 0] <- 0
pt.full$up[pt.full$up > 1] <- 1


test = group_by(dat,Start,StratifLength,IncubationTemp,dLabel) %>%
summarise(tGrm = sum(tGrm)/Seeds.Dish)
test$time = test$Start
test$pct = test$tGrm

#Plot of loglinear AFT, compare to KM curves
ggplot(test, aes(x=time, y = pct)) + theme_classic() + 
geom_line(data = pt.full,aes(color = IncubationTemp),size=1) + 
#geom_ribbon(data = pt.full,aes(ymin=pct-se,ymax=pct+se,group = IncubationTemp),alpha=0.15)+ # I actually don't think these SEs are for the pct germinated but an error of the time to that quantile...
facet_wrap(aes(group = StratifLength)) + 
xlim(c(0,5))+
ylim(c(0,1))+
labs(x = "Time (Weeks)", y = "Germination Proportion")+
theme(panel.grid = element_blank(), panel.border = element_blank(), text = element_text(size = 18), legend.text=element_text(size=16),legend.title=element_text(size=16),
	axis.line=element_line(size=1),axis.text.y=element_text(size=16),axis.text.x=element_text(size=16),axis.title=element_text(size=18),
	plot.title = element_text(size=16), legend.position = "top")+
guides(color=guide_legend(title=expression("Incubation Temperature " (degree~C))))
ggsave("Marmoh_Exp1_T2E_LINESONLY.tiff",units="in",width=10,height=6,dpi=300)


#Diagnostic and model predictions show poor fit for parametric AFT models.
#Use NP method for interval censored data using "interval" package. 
#The following is a loop to perform pairwise logrank tests using ictest(). 
#P-values are corrected using the bonferroni correction for 30 tests. 

#WARNING!!! When using permutational methods like wsr.mc (within-subjects resampling Monte Carlo)
#this loop will take a long time. I've tried to fix that by reprogramming it to only do the test for
#the lower matrix which should cut the time in half. Doing the full matrix of comparisons, minus the diagonal,
#takes hours when using 9,999 permutations, and I just don't see any reason to use fewer than that. 

expt1.ic.km = icfit(so ~ StratifLength + IncubationTemp, data = dat.expand2.2)

dat.expand2.2$trt = with(dat.expand2.2,interaction(StratifLength,IncubationTemp))

pmat = matrix(0,nrow=6,ncol=6)
scores = matrix(NA,nrow=6,ncol=6)
sig = paste("Significant p-value with bonferroni correction =", 0.05/((6*5)/2))
res = list(pmat,sig)
rownames(pmat) = with(dat.expand2.2,interaction(expand.grid(levels(StratifLength),levels(IncubationTemp))))
colnames(pmat) = with(dat.expand2.2,interaction(expand.grid(levels(StratifLength),levels(IncubationTemp))))
rownames(scores) = with(dat.expand2.2,interaction(expand.grid(levels(StratifLength),levels(IncubationTemp))))
colnames(scores) = with(dat.expand2.2,interaction(expand.grid(levels(StratifLength),levels(IncubationTemp))))

for (i in 1:6){
for (j in i:6){
dat.g1 = dat.expand2.2[dat.expand2.2$trt == levels(dat.expand2.2$trt)[i],]
dat.g2 = dat.expand2.2[dat.expand2.2$trt == levels(dat.expand2.2$trt)[j],]
dat.g1g2 = rbind(dat.g1,dat.g2)


so = with(dat.g1g2,Surv(time = Start, time2 = End, type ="interval2"))

if (i == j){

pmat[i,j] = 1
scores[i,j] = NA

}else if (i != j){

test = ictest(so ~ trt, data = dat.g1g2,scores="wmw",method = "wsr.HLY",mcontrol = mControl(np=1e5-1,nmc=1e5-1))

pmat[i,j] = round(test$p.value,4)
pmat[j,i] = round(test$p.value,4)
scores[i,j] = test$U[1]
scores[j,i] = -test$U[1]
}
}}

res = list(test.scores = scores, p.values = pmat,Significant = sig)
capture.output(res,file="pairwise_logrank.txt")


########################################################################################
##################################### Analysis 2 #######################################
########################################################################################

#Read in data for experiment 2
dat2 = read.csv("Marmoh_Analysis 2_cumulative.csv",header=T)

#Binomial GLM to examine cumulative germination
exp2.binom = glm(PropGerm ~ Temp2 * Temp2Photo, data = dat2, family = "binomial", weights = Seeds.Dish) #Main effects sig, but not interaction

#There is some overdispersion. Quasibinomial does not change estimates or results, but use because overdispersion is present
exp2.binom = glm(PropGerm ~ Temp2 * Temp2Photo, data = dat2, family = "quasibinomial", weights = Seeds.Dish)
capture.output(summary(exp2.binom), file = "Exp2_cumulative_GLM.txt")
exp2.binom.aov = Anova(exp2.binom, type = 3, test = "LR")
capture.output(exp2.binom.aov, file = "Exp2_cumulative_GLM_ANOVA.txt")

#Estimates of marginal means
exp2.binom.lsm = emmeans(exp2.binom, pairwise ~ Temp2 * Temp2Photo,type="response")

#Predicted means with group memberships
exp2.binom.pred = CLD(exp2.binom.lsm, Letters = LETTERS)
capture.output(exp2.binom.pred, file = "Exp2_cumulative_pred.txt")

#Plot results
ggplot(exp2.binom.pred, aes(x = Temp2, y = prob, group = Temp2Photo)) + theme_classic()+
geom_col(aes(fill = Temp2Photo),color = "Black",position = "dodge",lwd = 1.3)+
geom_errorbar(aes(ymin = prob - SE, ymax = prob + SE),position = position_dodge(0.9),width = 0.2, lwd = 1.3)+
scale_fill_manual(values = c("DarkGray","White"))+
labs(x = expression("Incubation Temperature " (degree~C)), y = "Germination Proportion")+
scale_y_continuous(expand=c(0,0),limits=c(0,1))+
theme(panel.grid = element_blank(), panel.border = element_blank(), text = element_text(size = 18), legend.text=element_text(size=16),legend.title=element_text(size=16),
	axis.line=element_line(size=1),axis.text.y=element_text(size=16),axis.text.x=element_text(size=16),axis.title=element_text(size=18),
	plot.title = element_text(size=16), legend.position = "top")+
guides(fill=guide_legend(title="Incubation Photoperiod"))
ggsave("Marmoh_Exp2_cumulative.tiff",units="in",width=6,height=6,dpi=600)

########################################################################################
##################################### Analysis 3 #######################################
########################################################################################
#Read data
dat3 = read.csv("Marmoh_Analysis 3_cumulative.csv",header=T)

#Change year to a factor, and set levels
dat3$Year = as.factor(dat3$Year)
levels(dat3$Year) <- c("2014","2017")

#Binomial GLM for cumulative germination
exp3.binom = glm(PropGerm ~ StartTemp * as.factor(Year), data = dat3, family = "binomial", weights = Seeds.Dish)	

#Slight overdispersion, quasibinomial does not change estimates or results, but use because overdispersion is present
exp3.binom = glm(PropGerm ~ StartTemp * as.factor(Year), data = dat3, family = "quasibinomial", weights = Seeds.Dish)
capture.output(summary(exp3.binom), file = "Exp3_cumulative_GLM.txt")
exp3.binom.aov = Anova(exp3.binom, type = 3, test="LR") #Apparently there is an effect by year..
capture.output(exp3.binom.aov, file = "Exp3_cumulative_GLM_ANOVA.txt")

#Estimate marginal means 
exp3.binom.lsm = emmeans(exp3.binom, pairwise ~ StartTemp * as.factor(Year), type="response")

#Predicted means with group memberships
exp3.binom.pred = CLD(exp3.binom.lsm,Letters = LETTERS)
capture.output(exp3.binom.pred, file = "Exp3_cumulative_pred.txt")

#Make sure that Year is a factor
exp3.binom.pred$Year = as.factor(exp3.binom.pred$Year)

#Plot results
ggplot(exp3.binom.pred, aes(x = StartTemp, y = prob, group = Year)) + theme_classic()+
geom_col(aes(fill = Year),color = "Black",position = "dodge",lwd = 1.3)+
geom_errorbar(aes(ymin = prob - SE, ymax = prob + SE),position = position_dodge(0.9),width = 0.2, lwd = 1.3)+
scale_fill_manual(values = c("DarkGray","White"))+
labs(x = expression("Incubation Temperature " (degree~C)), y = "Germination Proportion")+
scale_y_continuous(expand=c(0,0),limits=c(0,1))+
theme(panel.grid = element_blank(), panel.border = element_blank(), text = element_text(size = 18), legend.text=element_text(size=16),legend.title=element_text(size=16),
	axis.line=element_line(size=1),axis.text.y=element_text(size=16),axis.text.x=element_text(size=16),axis.title=element_text(size=18),
	plot.title = element_text(size=16), legend.position = "top")+
guides(fill=guide_legend(title="Year"))
ggsave("Marmoh_Exp3_cumulative.tiff",units="in",width=6,height=6,dpi=600)

########################################################################################
##################################### Analysis 4 #######################################
########################################################################################

#Read data
dat4 = read.csv("Marmoh_Analysis 4_cumulative.csv",header=T)

#Because there is total separation (no germination at sub-optimal temps in the dark) we have to use bias-reduced glm (brglmFit) from the brglm2 package. 
exp4.binom = glm(PropGerm ~ StartTemp * StartPhotoperiod, data = dat4, family = "binomial", weights = Seeds.Dish,method="brglmFit") 
capture.output(summary(exp4.binom), file = "Exp4_cumulative_GLM.txt")
exp4.binom.aov = Anova(exp4.binom, type = 3, test = "LR")
capture.output(exp4.binom.aov, file = "Exp4_cumulative_GLM_ANOVA.txt")

#Estimate marginal means
exp4.binom.lsm = emmeans(exp4.binom, pairwise ~ StartTemp * StartPhotoperiod, type = "response")

#Predicted means with group memberships
exp4.binom.pred = CLD(exp4.binom.lsm, Letters = LETTERS)
capture.output(exp4.binom.pred, file = "Exp4_cumulative_pred.txt")

#Plot results
ggplot(exp4.binom.pred, aes(x = StartTemp, y = prob, group = StartPhotoperiod)) + theme_classic()+
geom_col(aes(fill = StartPhotoperiod),color = "Black",position = "dodge",lwd = 1.3)+
geom_errorbar(aes(ymin = prob, ymax = prob + SE),position = position_dodge(0.9),width = 0, lwd = 1.3)+
scale_fill_manual(values = c("DarkGray","White"))+
labs(x = expression("Incubation Temperature " (degree~C)), y = "Percent Germinated")+
scale_y_continuous(expand=c(0,0),limits=c(0,1))+
theme(panel.grid = element_blank(), panel.border = element_blank(), text = element_text(size = 18), legend.text=element_text(size=16),legend.title=element_text(size=16),
	axis.line=element_line(size=1),axis.text.y=element_text(size=16),axis.text.x=element_text(size=16),axis.title=element_text(size=18),
	plot.title = element_text(size=16), legend.position = "top")+
guides(fill=guide_legend(title="Incubation Photoperiod"))
ggsave("Marmoh_Exp4_cumulative.tiff",units="in",width=6,height=6,dpi=600)
