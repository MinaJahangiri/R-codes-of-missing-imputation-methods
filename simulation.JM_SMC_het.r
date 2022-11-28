start.simulation<-Sys.time()

##########################
#Load packages
##########################
set.seed(1234)
library(purrr)
library(lme4)
library(REEMtree)
library(longitudinalData)
library(plyr)
library(jomo)
library(mitml)
library(mitools)
library(MuMIn)
library(merTools)
library(foreign)
M<-read.csv("simulation.keith.long.1000.csv",header=TRUE)
nrow(M)
names(M)
names(M)[names(M) == "period"] <- "time"
n.sim<-1000

##########################
#Imputation 
##########################
A<-NULL
for (i in 1:n.sim){ 
formula<-as.formula(dbp~age+sex+bmi+(1|id))
B<-jomo.smc(formula,data=M[((i-1)*nrow(M[which(M$dataset==i),])+1):((i-1)*nrow(M[which(M$dataset==i),])+nrow(M[which(M$dataset==i),])),]
, nburn=1000, nbetween=1000, nimp=10, meth="random",family="gaussian",model="lmer")
A<-rbind(B,A)
}
write.table(A,"simulation.jomo.smc.het.csv",sep=",",col.names=T,row.names=F)
A<-read.csv("simulation.jomo.smc.het.csv",header=TRUE, na.strings ="NA")
names(A)
A<-A[,-8]
names(A)[names(A) == "clus"] <- "id"
names(A)[names(A) == "Imputation"] <- "imp"
names(A)
str(A)
A<-A[which(A$imp!=0),]
A<-within(A, sex<-factor(sex))  # sex is a factor

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
  result1_df[i,]<-fixef(mod,add.dropped=TRUE) #extract coefficients to dataframe
rownames(result1_df)[i]<-names(list_df)[i] #assign rowname to results from data used
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
rownames(result2_df)[i]<-names(list_df)[i] #assign rowname to results from data used
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
rownames(result3_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result3_df
mean(result3_df$MSE)

###################RMSE for each imputation of lmer
result4_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result4_df)<-c("RMSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age + sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
    result4_df[i,]<- sqrt(mean((residuals(mod,type="pearson",scaled=TRUE))^2)) #extract RMSE to dataframe
rownames(result4_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result4_df
mean(result4_df$RMSE)

###################MAD for each imputation of lmer
result5_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result5_df)<-c("MAD") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age + sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
    result5_df[i,]<- mean(abs(residuals(mod,type="pearson",scaled=TRUE))) #extract MAD to dataframe
rownames(result5_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result5_df
mean(result5_df$MAD)

###################AIC for each imputation of lmer
result6_df <- as.data.frame(matrix(ncol=2,nrow=length(list_df))) # make an empty datafram
colnames(result6_df)<-c("df","AIC") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age +sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
    result6_df[i,]<- extractAIC(mod) #extract AIC to dataframe
rownames(result6_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
matrix(c(result6_df[,2]))
mean(result6_df[,2])

###################Deviance for each imputation of lmer
result7_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result7_df)<-c("Deviance") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age +sex+ bmi + time+(1|id),data=list_df[[i]]) #mixed model
    result7_df[i,]<- -2*((summary(mod))$logLik) #extract Deviance to dataframe
rownames(result7_df)[i]<-names(list_df)[i] #assign rowname to results from data used
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
rownames(result8_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result8_df
mean(result8_df$MSE)

###################RMSE for each imputation of REEMtrees
result9_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result9_df)<-c("RMSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
mod<-REEMtree(dbp~age +sex+ bmi + time,random=~1|id,data=na.omit(list_df[[i]]),lme.control=lmeControl(opt ="optim"))
 result9_df[i,]<- sqrt(mean((residuals(mod))^2)) #extract RMSE to dataframe
rownames(result9_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result9_df
mean(result9_df$RMSE)

###################MAD for each imputation of REEMtrees
result10_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result10_df)<-c("MAD") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
mod<-REEMtree(dbp~age +sex+ bmi + time,random=~1|id,data=na.omit(list_df[[i]]),lme.control=lmeControl(opt ="optim"))
    result10_df[i,]<- mean(abs(residuals(mod))) #extract MAD to dataframe
rownames(result10_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result10_df
mean(result10_df$MAD)

###################Deviance for each imputation of REEMtrees
result11_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result11_df)<-c("Deviance") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
mod<-REEMtree(dbp~age +sex+ bmi + time,random=~1|id,data=na.omit(list_df[[i]]),lme.control=lmeControl(opt ="optim"))
    result11_df[i,]<- -2*(logLik.REEMtree(mod)) #extract Deviance to dataframe
rownames(result11_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result11_df
mean(result11_df$Deviance)

end.simulation<-Sys.time()
time.simulation<- end.simulation-start.simulation
time.simulation



