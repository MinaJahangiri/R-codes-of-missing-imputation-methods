start.simulation<-Sys.time()

set.seed(1234)

##########################
#Load packages and data
#########################

library(lme4)
library(REEMtree)
library(MuMIn)
library(merTools)
library(longitudinalData)
M<-read.csv("simulation.keith.wide.1000.csv",header=TRUE)
n.sim<-1000

##########################
#split dataset 
##########################

list_df <- split(M, M$dataset) 

##########################
#single imputation
##########################

#bmi
copy.mean2.bmi<-NULL
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
bmi.traj<-as.matrix(M[((i-1)*nrow(M[which(M$dataset==i),])+1):((i-1)*nrow(M[which(M$dataset==i),])+nrow(M[which(M$dataset==i),])),c(4,8,11,14,17,20)])
imp.bmi<-imputation(bmi.traj,method="crossHotDeck",lowerBound="globalMin",upperBound="globalMax")
copy.mean1.bmi<-as.data.frame(imp.bmi)
copy.mean2.bmi<-rbind(copy.mean2.bmi,copy.mean1.bmi)
}

#dbp
copy.mean2.dbp<-NULL
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
dbp.traj<-as.matrix(M[((i-1)*nrow(M[which(M$dataset==i),])+1):((i-1)*nrow(M[which(M$dataset==i),])+nrow(M[which(M$dataset==i),])),c(5,9,12,15,18,21)])
imp.dbp<-imputation(dbp.traj,method="crossHotDeck",lowerBound="globalMin",upperBound="globalMax")
copy.mean1.dbp<-as.data.frame(imp.dbp)
copy.mean2.dbp<-rbind(copy.mean2.dbp,copy.mean1.dbp)
}

###################long data
data<-cbind(M[,c(1:3,6,7,10,13,16,19)],copy.mean2.bmi,copy.mean2.dbp)
nrow(data)
names(data)
summary(data)

A<-NULL
for (i in 1:n.sim){
B<-reshape(data[((i-1)*nrow(data[which(data$dataset==i),])+1):((i-1)*nrow(data[which(data$dataset==i),])+nrow(data[which(data$dataset==i),])),],varying=list(c("age0","age1","age2","age3","age4","age5"),c("bmi0","bmi1","bmi2","bmi3","bmi4","bmi5"),
c("dbp0","dbp1","dbp2","dbp3","dbp4","dbp5")),idvar="data$id",
v.names=c("age","bmi","dbp"),times=c(0,1,2,3,4,5),direction="long")
A<-rbind(B,A)
}
names(A)
nrow(A)

##########################
#split dataset by imp
##########################

list_df <- split(A, A$dataset) 

################################################
#lmer
###############################################
###################Coeficient for each imputation of lmer
result1_df <- as.data.frame(matrix(ncol=5,nrow=length(list_df))) # make an empty dataframe
colnames(result1_df)<-c("intercept","age","sex","bmi","time") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age +sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
  result1_df[i,]<-fixef(mod) #extract coefficients to dataframe
}
result1_df
mean(result1_df$intercept)
mean(result1_df$age)
mean(result1_df$sex)
mean(result1_df$bmi)
mean(result1_df$time)

###################StD.Error for each imputation of lmer
result2_df <- as.data.frame(matrix(ncol=5,nrow=length(list_df))) # make an empty datafram
colnames(result2_df)<-c("intercept","age","sex","bmi","time") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age + sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
    result2_df[i,]<- se.fixef(mod) #extract StD.Error to dataframe
}
result2_df
mean(result2_df$intercept)
mean(result2_df$age)
mean(result2_df$sex)
mean(result2_df$bmi)
mean(result2_df$time)

###################MSE for each imputation of lmer
result3_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result3_df)<-c("MSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age +sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
    result3_df[i,]<- mean((residuals(mod,type="pearson",scaled=TRUE))^2) #extract MSE to dataframe
}
result3_df
mean(result3_df$MSE)

###################RMSE for each imputation of lmer
result4_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result4_df)<-c("RMSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age + sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
    result4_df[i,]<- sqrt(mean((residuals(mod,type="pearson",scaled=TRUE))^2)) #extract RMSE to dataframe
}
result4_df
mean(result4_df$RMSE)

###################MAD for each imputation of lmer
result5_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result5_df)<-c("MAD") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age + sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
    result5_df[i,]<- mean(abs(residuals(mod,type="pearson",scaled=TRUE))) #extract MAD to dataframe
}
result5_df
mean(result5_df$MAD)

###################AIC for each imputation of lmer
result6_df <- as.data.frame(matrix(ncol=2,nrow=length(list_df))) # make an empty datafram
colnames(result6_df)<-c("df","AIC") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age +sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
    result6_df[i,]<- extractAIC(mod) #extract AIC to dataframe
}
matrix(c(result6_df[,2]))
mean(result6_df[,2])

###################Deviance for each imputation of lmer
result7_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result7_df)<-c("Deviance") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age +sex+ bmi + time+(1|id),data=list_df[[i]]) #mixed model
    result7_df[i,]<- -2*((summary(mod))$logLik) #extract Deviance to dataframe
}
result7_df
mean(result7_df$Deviance)

################################################
#REEM trees
###############################################
###################MSE for each imputation of REEMtrees
result8_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result8_df)<-c("MSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
mod<-REEMtree(dbp~age +sex+ bmi + time,random=~1|id,data=na.omit(list_df[[i]]),lme.control=lmeControl(opt ="optim"))
 result8_df[i,]<- mean((residuals(mod))^2) #extract MSE to dataframe
}
result8_df
mean(result8_df$MSE)

###################RMSE for each imputation of REEMtrees
result9_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result9_df)<-c("RMSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
mod<-REEMtree(dbp~age +sex+ bmi + time,random=~1|id,data=na.omit(list_df[[i]]),lme.control=lmeControl(opt ="optim"))
 result9_df[i,]<- sqrt(mean((residuals(mod))^2)) #extract RMSE to dataframe
}
result9_df
mean(result9_df$RMSE)

###################MAD for each imputation of REEMtrees
result10_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result10_df)<-c("MAD") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
mod<-REEMtree(dbp~age +sex+ bmi + time,random=~1|id,data=na.omit(list_df[[i]]),lme.control=lmeControl(opt ="optim"))
    result10_df[i,]<- mean(abs(residuals(mod))) #extract MAD to dataframe
}
result10_df
mean(result10_df$MAD)

###################Deviance for each imputation of REEMtrees
result11_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result11_df)<-c("Deviance") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
mod<-REEMtree(dbp~age +sex+ bmi + time,random=~1|id,data=na.omit(list_df[[i]]),lme.control=lmeControl(opt ="optim"))
    result11_df[i,]<- -2*(logLik.REEMtree(mod)) #extract Deviance to dataframe
}
result11_df
mean(result11_df$Deviance)

end.simulation<-Sys.time()
time.simulation<- end.simulation-start.simulation
time.simulation


