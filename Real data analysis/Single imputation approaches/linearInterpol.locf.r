rm(list=ls())
#########Load packages
library(lme4)
library(geepack)
library(REEMtree)
library(misty)
library(MuMIn)
library(merTools)
library(longitudinalData)
set.seed(1234)

#######################################
# Data preparation 
#######################################
setwd("C:/Users/jahangiri/Desktop")
A<-read.csv("A.csv",header=TRUE)
A[A==-9]<-NA
names(A)
B<-A[,c(1,22:28,57:62,108:113,173:178)]
names(B)
n<-length(B$id);n
str(B)
B<-within(B, sex<-factor(Gender))  # gender is a factor
B$id<-as.numeric(as.factor(B$id))
str(B)
B<-B[order(B$id),]
head(B)
tail(B)

##########################################################
# Imputation for trajectory of bmi variable using linearInterpol.locf
##########################################################
bmi<-B[,3:8]
names(bmi)
head(bmi)
tail(bmi)
bmi.traj<-as.matrix(bmi)
imp.bmi<-imputation(bmi.traj,method="linearInterpol.locf",lowerBound="globalMin",upperBound="globalMax")
bmi<-as.data.frame(imp.bmi)
colnames(bmi)<-c("bmi1","bmi2","bmi3","bmi4","bmi5","bmi6")
head(bmi)
tail(bmi)

##########################################################
# Imputation for trajectory of dbp variable using linearInterpol.locf
##########################################################
dbp<-B[,9:14]
names(dbp)
head(dbp)
tail(dbp)
dbp.traj<-as.matrix(dbp)
imp.dbp<-imputation(dbp.traj,method="linearInterpol.locf",lowerBound="globalMin",upperBound="globalMax")
dbp<-as.data.frame(imp.dbp)
colnames(dbp)<-c("dbp1","dbp2","dbp3","dbp4","dbp5","dbp6")
head(dbp)
tail(dbp)

##########################################################
# Imputation for trajectory of sbp variable using linearInterpol.locf
##########################################################
sbp<-B[,15:20]
names(sbp)
head(sbp)
tail(sbp)
sbp.traj<-as.matrix(sbp)
imp.sbp<-imputation(sbp.traj,method="linearInterpol.locf",lowerBound="globalMin",upperBound="globalMax")
sbp<-as.data.frame(imp.sbp)
colnames(sbp)<-c("sbp1","sbp2","sbp3","sbp4","sbp5","sbp6")
head(sbp)
tail(sbp)

#######################################################
#cbind imputed trajectories and transform to long format
######################################################
age<-B[,21:26]
names(age)
B<-cbind(B$id,B$sex,age,bmi,dbp,sbp)
names(B)
names(B)[names(B) == "B$id"] <- "id"
names(B)[names(B) == "B$sex"] <- "sex"
names(B)
str(B)

#reshape to long
A<-reshape(B,varying=list(c("age1","age2","age3","age4","age5","age6"),c("bmi1","bmi2","bmi3","bmi4","bmi5","bmi6")
,c("dbp1","dbp2","dbp3","dbp4","dbp5","dbp6"),c("sbp1","sbp2","sbp3","sbp4","sbp5","sbp6")),idvar="id",
v.names=c("age","bmi","dbp","sbp"),times=c(1,2,3,4,5,6),direction="long")
names(A)
str(A)
head(A)
tail(A)

###########################################################
#Data analysis for DBP
#############################################################

#lmer
fit.lmer<-lmer(dbp ~ age+sex+bmi+time+(1|id),data=A)
summary(fit.lmer)

#multicollinearity
VIF<-collin.diag(fit.lmer, print = c("all", "vif", "eigen"), digits = 3, p.digits = 3,
check = TRUE, output = TRUE);VIF

MSE<- mean((residuals(fit.lmer,type="pearson",scaled=TRUE))^2);MSE 
RMSE<- sqrt(mean((residuals(fit.lmer,type="pearson",scaled=TRUE))^2));RMSE
MAD<- mean(abs(residuals(fit.lmer,type="pearson",scaled=TRUE)));MAD
AIC<- extractAIC(fit.lmer);AIC
Deviance<- -2*((summary(fit.lmer))$logLik);Deviance

#tree
fit.tree<-REEMtree(dbp ~ age+sex+ bmi,random=~1|id,data=A,lme.control=lmeControl(opt ="optim"))
residuals.tree<-residuals(fit.tree)
MSE<-mean((residuals.tree)^2);MSE
RMSE<-sqrt(mean((residuals.tree)^2));RMSE
MAD<- mean(abs(residuals.tree));MAD
Deviance<- -2*(logLik.REEMtree(fit.tree));Deviance
plot.REEMtree(fit.tree,text=TRUE)
tree(fit.tree)

###########################################################
#Data analysis for SBP
#############################################################

#lmer
fit.lmer<-lmer(sbp ~ age+sex+bmi+time+(1|id),data=A)
summary(fit.lmer)

MSE<- mean((residuals(fit.lmer,type="pearson",scaled=TRUE))^2);MSE 
RMSE<- sqrt(mean((residuals(fit.lmer,type="pearson",scaled=TRUE))^2));RMSE
MAD<- mean(abs(residuals(fit.lmer,type="pearson",scaled=TRUE)));MAD
AIC<- extractAIC(fit.lmer);AIC
Deviance<- -2*((summary(fit.lmer))$logLik);Deviance

#tree
fit.tree<-REEMtree(sbp ~ age+sex+ bmi,random=~1|id,data=A,lme.control=lmeControl(opt ="optim"))
residuals.tree<-residuals(fit.tree)
MSE<-mean((residuals.tree)^2);MSE
RMSE<-sqrt(mean((residuals.tree)^2));RMSE
MAD<- mean(abs(residuals.tree));MAD
Deviance<- -2*(logLik.REEMtree(fit.tree));Deviance
plot.REEMtree(fit.tree,text=TRUE)
tree(fit.tree)


