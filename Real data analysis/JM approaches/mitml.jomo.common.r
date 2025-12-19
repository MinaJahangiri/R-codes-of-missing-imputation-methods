#########Load packages
library(jomo)
library(mitml)
library(mitools)
library(lme4)
library(geepack)
library(REEMtree)
library(MuMIn)
library(merTools)
library(akmedoids)
library(longitudinalData)
library(foreign)
set.seed(1234)

#######################################
# Data preparation
#######################################
A<-read.csv("A.csv",header=TRUE)
A[A==-9]<-NA
names(A)
A<-A[,c(1,22:28,57:62,173:178)]
names(A)

#reshape to long
A<-reshape(A,varying=list(c("age1","age2","age3","age4","age5","age6"),c("bmi1","bmi2","bmi3","bmi4","bmi5","bmi6"),
c("dbp1","dbp2","dbp3","dbp4","dbp5","dbp6")),idvar="id",
v.names=c("age","bmi","dbp"),times=c(1,2,3,4,5,6),direction="long")
names(A)
str(A)
head(A)
tail(A)
A<-within(A, sex<-factor(Gender))   # gender is a factor
A$id<-as.numeric(as.factor(A$id))
str(A)
A<-A[order(A$id),]

############################################
# Multiple imputation using mitml.jomo.common
############################################
formula<-dbp+age+bmi+sex~1+(1|id)

imp<-jomoImpute(data=A, formula=formula,random.L1="none", n.burn=5000,
n.iter=1000, m=20,seed=1234)

#convergency
summary(imp)

pdf("mitml.jomo.common.DBP.export.pdf")
plot(imp,trace="all",export="pdf")   # trace="imputation"
dev.off()

pdf("mitml.jomo.common.DBP.export.pdf")
plot(imp,trace="all",export="none")   # trace="imputation"
dev.off()

#######################################
#write data
#######################################
#a list containing the imputed data sets
implist <- mitmlComplete(imp, print="all")
#write data
A<-write.mitmlSAV(implist, filename="mitml.jomo.common.DBP.sav")
A<-write.table(read.spss("mitml.jomo.common.DBP.sav"), file="mitml.jomo.common.DBP.csv", quote = TRUE, col.names=T,row.names=F,sep = ",")
B<-read.csv("mitml.jomo.common.DBP.csv",header=TRUE, na.string ="NA")
names(B)
str(B)
B<-B[ , -7] 
B<-within(B, sex<-factor(sex))   # sex is a factor
str(B)
names(B)[names(B) == "Imputation_"] <- "imp"
B$imp<-B$imp+1
names(B)
head(B)
tail(B)

##########################
#split dataset by imp
##########################
list_df <- split(B, B$imp) 

################################################
#lmer
###############################################
###################Coeficient for each imputation of lmer
result1_df <- as.data.frame(matrix(ncol=5,nrow=length(list_df))) # make an empty dataframe
colnames(result1_df)<-c("intercept","age","sex","bmi","time") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age +sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
  result1_df[i,]<-fixef(mod) #extract coefficients to dataframe
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
  mod<-lmer(dbp~age +sex+ bmi + time+(1|id),data=list_df[[i]]) #mixed model
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
  mod<-lmer(dbp~age +sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
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
  mod<-lmer(dbp~age +sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
    result7_df[i,]<- -2*((summary(mod))$logLik) #extract Deviance to dataframe
rownames(result7_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result7_df
mean(result7_df$Deviance)

################################################
#REEM trees
###############################################
###################MSE for each imputation of REEMtrees
result3_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result3_df)<-c("MSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
mod<-REEMtree(dbp~age +sex+ bmi ,random=~1|id,data=list_df[[i]],lme.control=lmeControl(opt ="optim"))
 result3_df[i,]<- mean((residuals(mod))^2) #extract MSE to dataframe
rownames(result3_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result3_df
mean(result3_df$MSE)
which(result3_df==min(result3_df))

###################RMSE for each imputation of REEMtrees
result4_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result4_df)<-c("RMSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
mod<-REEMtree(dbp~age +sex+ bmi ,random=~1|id,data=list_df[[i]],lme.control=lmeControl(opt ="optim"))
 result4_df[i,]<- sqrt(mean((residuals(mod))^2)) #extract RMSE to dataframe
rownames(result4_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result4_df
mean(result4_df$RMSE)
which(result4_df==min(result4_df))

###################MAD for each imputation of REEMtrees
result5_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result5_df)<-c("MAD") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
mod<-REEMtree(dbp~age +sex+ bmi ,random=~1|id,data=list_df[[i]],lme.control=lmeControl(opt ="optim"))
    result5_df[i,]<- mean(abs(residuals(mod))) #extract MAD to dataframe
rownames(result5_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result5_df
mean(result5_df$MAD)
which(result5_df==min(result5_df))

###################Deviance for each imputation of REEMtrees
result6_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result6_df)<-c("Deviance") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
mod<-REEMtree(dbp~age +sex+ bmi ,random=~1|id,data=list_df[[i]],lme.control=lmeControl(opt ="optim"))
    result6_df[i,]<- -2*(logLik.REEMtree(mod)) #extract Deviance to dataframe
rownames(result6_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result6_df

mean(result6_df$Deviance)
which(result6_df==min(result6_df))

