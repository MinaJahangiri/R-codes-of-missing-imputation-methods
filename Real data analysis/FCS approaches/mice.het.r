#########Load packages
library(lme4)
library(geepack)
library(REEMtree)
library(MuMIn)
library(merTools)
library(mice)
library(micemd)
library(miceadds)
library(mitools)
library(broom.mixed)
library(akmedoids)
library(longitudinalData)

#######################################
# Data preparation
#######################################
set.seed(1234)
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
A<-A[,-2]
names(A)

################################################
# empty imputation in mice
################################################
imp0 <- mice(A, maxit=0)
predM <- imp0$predictorMatrix
impM <- imp0$method

###############################################
# specify predictor matrix and imputation Method
##############################################
predM1 <- predM
predM1["dbp",] <- c(-2,0,1,1,0,1)
predM1["time",] <- c(-2,0,0,0,0,0)
predM1["id",] <- c(0,0,0,0,0,0)
predM1["age",] <- c(-2,0,0,0,0,0)
predM1["bmi",] <- c(-2,0,1,0,1,1)
predM1["sex",] <- c(-2,0,0,0,0,0)
predM1

impM1 <- impM
impM1["dbp"] <- "2l.norm"    #2l.lmer , #2l.norm (heter) , #CART and RF for other variables
impM1["bmi"] <- "2l.norm"
impM1

###########################################
# Multiple imputation using package mice
###########################################
imp <- mice(A, m = 20, predictorMatrix = predM1,Method = impM1, maxit=20,diagnostics=TRUE,seed=1234)

#######################################
# plots
#######################################
pdf("mice.het.impplot.DBP.pdf")
plot(imp)
dev.off()

pdf("mice.het.bwplot.DBP.pdf")
bwplot(imp)
dev.off()

pdf("mice.het.stripplot.DBP.pdf")
stripplot(imp)
dev.off()

pdf("mice.het.densityplot.DBP.pdf")
densityplot(imp)
dev.off()

pdf("mice.het.densityplotdbp.DBP.pdf")
densityplot(imp, ~dbp|.imp)   #blue is observed, red is imputed
dev.off()

pdf("mice.het.densityplotbmi.DBP.pdf")
densityplot(imp, ~bmi|.imp)   #blue is observed, red is imputed
dev.off()

#######################################
#Complete data and import to csv file
#######################################
data<-complete(imp, action = "long", include = TRUE)
write.table(data,"mice.het.DBP.csv",sep=",",col.names=T,row.names=F)
B<-read.csv("mice.het.DBP.csv",header=TRUE,na.string="NA")
names(B)
str(B)
B<-within(B, sex<-factor(sex))  # sex is a factor
B$age<-as.numeric(B$age)
str(B)
names(B)[names(B) == ".imp"] <- "imp"
B<-B[ , -2] 
names(B)
B<-B[which(B$imp!=0),]
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

######################################
#GEE
######################################
###################Coeficient for each imputation of GEE
result1_df <- as.data.frame(matrix(ncol=5,nrow=length(list_df))) # make an empty dataframe
colnames(result1_df)<-c("intercept","age","sex","bmi","time") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-geeglm(dbp~age + sex + bmi + time,id=id,family="gaussian",corstr="independence",data=list_df[[i]]) #GEE model
  result1_df[i,]<-mod$coefficients #extract coefficients to dataframe
rownames(result1_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result1_df
mean(result1_df$intercept)
mean(result1_df$age)
mean(result1_df$sex)
mean(result1_df$bmi)
mean(result1_df$time)

###################StD.Error for each imputation of GEE
result2_df <- as.data.frame(matrix(ncol=5,nrow=length(list_df))) # make an empty datafram
colnames(result2_df)<-c("intercept","age","sex","bmi","time") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-geeglm(dbp~age + sex + bmi + time,id=id,family="gaussian",corstr="independence",data=list_df[[i]]) #GEE model
    result2_df[i,]<- summary(mod)$coefficients[,2] #extract StD.Error to dataframe
rownames(result2_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result2_df
mean(result2_df$intercept)
mean(result2_df$age)
mean(result2_df$sex)
mean(result2_df$bmi)
mean(result2_df$time)

###################MSE for each imputation of GEE
result3_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result3_df)<-c("MSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-geeglm(dbp~age+sex+ bmi + time,id=id,family="gaussian",corstr="independence",data=list_df[[i]]) #GEE model
    result3_df[i,]<- mean((residuals(mod))^2) #extract MSE to dataframe
rownames(result3_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result3_df
mean(result3_df$MSE)

###################RMSE for each imputation of GEE
result4_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result4_df)<-c("RMSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-geeglm(dbp~age + sex + bmi + time,id=id,family="gaussian",corstr="independence",data=list_df[[i]]) #GEE model
    result4_df[i,]<- sqrt(mean((residuals(mod))^2)) #extract MSE to dataframe
rownames(result4_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result4_df
mean(result4_df$RMSE)

###################MAD for each imputation of GEE
result5_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result5_df)<-c("MAD") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-geeglm(dbp~age +sex + bmi + time,id=id,family="gaussian",corstr="independence",data=list_df[[i]]) #GEE model
    result5_df[i,]<- mean(abs(residuals(mod))) #extract MAD to dataframe
rownames(result5_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result5_df
mean(result5_df$MAD)

###################QIC for each imputation of GEE
result6_df <- as.data.frame(matrix(ncol=2,nrow=length(list_df))) # make an empty datafram
colnames(result6_df)<-c("df","AIC") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-geeglm(dbp~age +sex + bmi + time,id=id,family="gaussian",corstr="independence",data=list_df[[i]]) #GEE model
    result6_df[i,]<- QIC(mod,typeR=FALSE) #extract QIC to dataframe
rownames(result6_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
matrix(c(result6_df[,2]))

mean(result6_df[,2])

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

#############reshape to wide
A<-reshape(B,idvar=c("id","sex","imp"),timevar="time",direction="wide")
write.table(A,"wide.mice.het.DBP.csv",sep=",",col.names=T,row.names=F)
A<-read.csv("wide.mice.het.DBP.csv",header=TRUE)
names(A)
head(A)
tail(A)

###################################
#age
###################################
age<-A[age<-substr(names(A),1,3)%in%c("age")]
age<-cbind(age,A$imp)
names(age)
head(age)
tail(age)
names(age)[names(age) == "A$imp"] <- "imp"
list_df <- split(age,age$imp) 
age_df <- as.data.frame(matrix(ncol=6,nrow=length(A$id)/20)) # make an empty dataframe
colnames(age_df)<-c("age1","age2","age3","age4","age5","age6") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-outlier_detect(as.matrix(age[which(age$imp==i),]), method = 1, threshold = 0.95,count = 1, replace_with = 1)
age_df[[i]]<-mod$Outliers_Replaced #extract coefficients to dataframe
}

age<- do.call("rbind",list(
age_df[[1]][,-7], age_df[[2]][,-7],
age_df[[3]][,-7], age_df[[4]][,-7],
age_df[[5]][,-7], age_df[[6]][,-7],
age_df[[7]][,-7], age_df[[8]][,-7],
age_df[[9]][,-7], age_df[[10]][,-7],
age_df[[11]][,-7], age_df[[12]][,-7],
age_df[[13]][,-7], age_df[[14]][,-7],
age_df[[15]][,-7], age_df[[16]][,-7],
age_df[[17]][,-7], age_df[[18]][,-7],
age_df[[19]][,-7], age_df[[20]][,-7]))
head(age)
tail(age)

agelong<-reshapeWideToLong(cbind(rep(1:length(A$id)/20,20),age))
names(agelong)
names(agelong)[names(agelong) == "values"] <- "age"
age<-agelong$age
head(age)
tail(age)

###################################
#dbp
###################################
dbp<-A[dbp<-substr(names(A),1,3)%in%c("dbp")]
dbp<-cbind(dbp,A$imp)
names(dbp)
head(dbp)
tail(dbp)
names(dbp)[names(dbp) == "A$imp"] <- "imp"
list_df <- split(dbp,dbp$imp) 
dbp_df <- as.data.frame(matrix(ncol=6,nrow=length(A$id)/20)) # make an empty dataframe
colnames(dbp_df)<-c("dbp1","dbp2","dbp3","dbp4","dbp5","dbp6") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-outlier_detect(as.matrix(dbp[which(dbp$imp==i),]), method = 1, threshold = 0.95,count = 1, replace_with = 1)
dbp_df[[i]]<-mod$Outliers_Replaced #extract coefficients to dataframe
}

dbp<- do.call("rbind",list(
dbp_df[[1]][,-7], dbp_df[[2]][,-7],
dbp_df[[3]][,-7], dbp_df[[4]][,-7],
dbp_df[[5]][,-7], dbp_df[[6]][,-7],
dbp_df[[7]][,-7], dbp_df[[8]][,-7],
dbp_df[[9]][,-7], dbp_df[[10]][,-7],
dbp_df[[11]][,-7], dbp_df[[12]][,-7],
dbp_df[[13]][,-7], dbp_df[[14]][,-7],
dbp_df[[15]][,-7], dbp_df[[16]][,-7],
dbp_df[[17]][,-7], dbp_df[[18]][,-7],
dbp_df[[19]][,-7], dbp_df[[20]][,-7]))
head(dbp)
tail(dbp)

dbplong<-reshapeWideToLong(cbind(rep(1:length(A$id)/20,20),dbp))
names(dbplong)
names(dbplong)[names(dbplong) == "values"] <- "dbp"
dbp<-dbplong$dbp
head(dbp)
tail(dbp)

###################################
#bmi
###################################
bmi<-A[bmi<-substr(names(A),1,3)%in%c("bmi")]
bmi<-cbind(bmi,A$imp)
names(bmi)
head(bmi)
tail(bmi)
names(bmi)[names(bmi) == "A$imp"] <- "imp"
list_df <- split(bmi,bmi$imp) 
bmi_df <- as.data.frame(matrix(ncol=6,nrow=length(A$id)/20)) # make an empty dataframe
colnames(bmi_df)<-c("bmi1","bmi2","bmi3","bmi4","bmi5","bmi6") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-outlier_detect(as.matrix(bmi[which(bmi$imp==i),]), method = 1, threshold = 0.95,count = 1, replace_with = 1)
bmi_df[[i]]<-mod$Outliers_Replaced #extract coefficients to dataframe
}

bmi<- do.call("rbind",list(
bmi_df[[1]][,-7], bmi_df[[2]][,-7],
bmi_df[[3]][,-7], bmi_df[[4]][,-7],
bmi_df[[5]][,-7], bmi_df[[6]][,-7],
bmi_df[[7]][,-7], bmi_df[[8]][,-7],
bmi_df[[9]][,-7], bmi_df[[10]][,-7],
bmi_df[[11]][,-7], bmi_df[[12]][,-7],
bmi_df[[13]][,-7], bmi_df[[14]][,-7],
bmi_df[[15]][,-7], bmi_df[[16]][,-7],
bmi_df[[17]][,-7], bmi_df[[18]][,-7],
bmi_df[[19]][,-7], bmi_df[[20]][,-7]))
head(bmi)
tail(bmi)
names(bmi)
bmilong<-reshapeWideToLong(cbind(rep(1:length(A$id)/20,20),bmi))
names(bmilong)
names(bmilong)[names(bmilong) == "values"] <- "bmi"
bmi<-bmilong$bmi
head(bmi)
tail(bmi)

########################################
#final data set for modelling
########################################
imp<-B$imp
time<-B$time
sex<-B$sex
id<-B$id
A<-as.data.frame(cbind(imp,id,time,sex,age,bmi,dbp))
names(A)
str(A)
A<-within(A, sex<-factor(sex))# sex is a factor
head(A)
tail(A)

##########################
#split dataset by imp
##########################
list_df <- split(A, A$imp) 

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
    result3_df[i,]<- mean((residuals(mod))^2) #extract MSE to dataframe
rownames(result3_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result3_df
mean(result3_df$MSE)

###################RMSE for each imputation of lmer
result4_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result4_df)<-c("RMSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age + sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
    result4_df[i,]<- sqrt(mean((residuals(mod))^2)) #extract RMSE to dataframe
rownames(result4_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result4_df
mean(result4_df$RMSE)

###################MAD for each imputation of lmer
result5_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result5_df)<-c("MAD") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-lmer(dbp~age + sex + bmi + time+(1|id),data=list_df[[i]]) #mixed model
    result5_df[i,]<- mean(abs(residuals(mod))) #extract MAD to dataframe
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

######################################
#GEE
######################################
###################Coeficient for each imputation of GEE
result1_df <- as.data.frame(matrix(ncol=5,nrow=length(list_df))) # make an empty dataframe
colnames(result1_df)<-c("intercept","age","sex","bmi","time") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-geeglm(dbp~age + sex + bmi + time,id=id,family="gaussian",corstr="independence",data=list_df[[i]]) #GEE model
  result1_df[i,]<-mod$coefficients #extract coefficients to dataframe
rownames(result1_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result1_df
mean(result1_df$intercept)
mean(result1_df$age)
mean(result1_df$sex)
mean(result1_df$bmi)
mean(result1_df$time)

###################StD.Error for each imputation of GEE
result2_df <- as.data.frame(matrix(ncol=5,nrow=length(list_df))) # make an empty datafram
colnames(result2_df)<-c("intercept","age","sex","bmi","time") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-geeglm(dbp~age + sex + bmi + time,id=id,family="gaussian",corstr="independence",data=list_df[[i]]) #GEE model
    result2_df[i,]<- summary(mod)$coefficients[,2] #extract StD.Error to dataframe
rownames(result2_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result2_df
mean(result2_df$intercept)
mean(result2_df$age)
mean(result2_df$sex)
mean(result2_df$bmi)
mean(result2_df$time)

###################MSE for each imputation of GEE
result3_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result3_df)<-c("MSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-geeglm(dbp~age +sex + bmi + time,id=id,family="gaussian",corstr="independence",data=list_df[[i]]) #GEE model
    result3_df[i,]<- mean((residuals(mod))^2) #extract MSE to dataframe
rownames(result3_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result3_df
mean(result3_df$MSE)

###################RMSE for each imputation of GEE
result4_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result4_df)<-c("RMSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-geeglm(dbp~age + sex + bmi + time,id=id,family="gaussian",corstr="independence",data=list_df[[i]]) #GEE model
    result4_df[i,]<- sqrt(mean((residuals(mod))^2)) #extract MSE to dataframe
rownames(result4_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result4_df
mean(result4_df$RMSE)

###################MAD for each imputation of GEE
result5_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result5_df)<-c("MAD") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-geeglm(dbp~age + sex + bmi + time,id=id,family="gaussian",corstr="independence",data=list_df[[i]]) #GEE model
    result5_df[i,]<- mean(abs(residuals(mod))) #extract MAD to dataframe
rownames(result5_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result5_df
mean(result5_df$MAD)

###################QIC for each imputation of GEE
result6_df <- as.data.frame(matrix(ncol=2,nrow=length(list_df))) # make an empty datafram
colnames(result6_df)<-c("df","AIC") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
  mod<-geeglm(dbp~age + sex + bmi + time,id=id,family="gaussian",corstr="independence",data=list_df[[i]]) #GEE model
    result6_df[i,]<- QIC(mod,typeR=FALSE) #extract QIC to dataframe
rownames(result6_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
matrix(c(result6_df[,2]))
mean(result6_df[,2])

################################################
#REEM trees
###############################################
###################MSE for each imputation of REEMtrees
result3_df <- as.data.frame(matrix(ncol=1,nrow=length(list_df))) # make an empty datafram
colnames(result3_df)<-c("MSE") #give the dataframe column names
for (i in 1:length(list_df)){ #run a loop over the dataframes in the list
mod<-REEMtree(dbp~age +sex+ bmi + time,random=~1|id,data=list_df[[i]],lme.control=lmeControl(opt ="optim"))
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
mod<-REEMtree(dbp~age +sex+ bmi + time,random=~1|id,data=list_df[[i]],lme.control=lmeControl(opt ="optim"))
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
mod<-REEMtree(dbp~age +sex+ bmi + time,random=~1|id,data=list_df[[i]],lme.control=lmeControl(opt ="optim"))
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
mod<-REEMtree(dbp~age +sex+ bmi + time,random=~1|id,data=list_df[[i]],lme.control=lmeControl(opt ="optim"))
    result6_df[i,]<- -2*(logLik.REEMtree(mod)) #extract Deviance to dataframe
rownames(result6_df)[i]<-names(list_df)[i] #assign rowname to results from data used
}
result6_df
mean(result6_df$Deviance)
which(result6_df==min(result6_df))





