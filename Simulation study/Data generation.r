start.simulation<-Sys.time()
set.seed(1234)
library(simstudy)
library(data.table)
library(msm)
n.sim<-1000
n<-1000
data<-NULL
for (i in 1:n.sim){

trunc_norm <- function(mean, sd, lower, upper) {
  rtnorm(n = 1, mean = mean, sd = sd, lower = lower, upper = upper)
}

defI <- defData(varname = "sex", 
  formula = 0.5, dist = "binary")
defI <- defData(defI, varname = "base_age", 
  formula = "trunc_norm(..mean, ..sd, ..lower, ..upper)",
  dist = "nonrandom")
defI <- defData(defI, varname = "theta_bmi", 
  formula = 0, variance = 4.50^2, dist = "normal")
defI <- defData(defI, varname = "theta_dbp", 
  formula = 0, variance = 6.35^2, dist = "normal")

defP <- defDataAdd(varname = "age", formula = "base_age + period*3",
  dist = "nonrandom")
defP <- defDataAdd(defP, varname = "bmi", 
  formula = "19.86 + 0.136*age + 2.360 * sex + theta_bmi",
  variance = 1.57 ^ 2, dist = "normal")
defP <- defDataAdd(defP, varname = "dbp", 
  formula = "55.03 + 0.098 * age - 3.434 * sex + 0.707 * bmi + theta_dbp",
  variance = 6.88^2, dist = "normal")

mean <- 39.34
sd <- 16.23
lower <- 1
upper <- 84

dd <- genData(n, defI)

dp <- addPeriods(dd, nPeriods = 6)
dp <- addColumns(defP, dp)

defB <- defCondition(condition = "period == 0", formula = -2.646, dist = "nonrandom")
defB <- defCondition(defB, "period == 1", formula = -2.634, dist = "nonrandom")
defB <- defCondition(defB, "period == 2", formula = -3.047, dist = "nonrandom")
defB <- defCondition(defB, "period == 3", formula = -3.271, dist = "nonrandom")
defB <- defCondition(defB, "period == 4", formula = -2.872, dist = "nonrandom")
defB <- defCondition(defB, "period == 5", formula = -2.440, dist = "nonrandom")

defD <- defCondition(condition = "period == 0", formula = -1.701, dist = "nonrandom")
defD <- defCondition(defD, "period == 1", formula = -1.815, dist = "nonrandom")
defD <- defCondition(defD, "period == 2", formula = -2.091, dist = "nonrandom")
defD <- defCondition(defD, "period == 3", formula = -2.386, dist = "nonrandom")
defD <- defCondition(defD, "period == 4", formula = -2.104, dist = "nonrandom")
defD <- defCondition(defD, "period == 5", formula = -1.522, dist = "nonrandom")

defM<- defMiss(varname = "bmi", 
  formula = "int_B + 0.002 * age + 0.02 * dbp", logit.link = TRUE)
defM<- defMiss(defM, varname = "dbp", 
  formula = "int_D + 0.002 * age + 0.02 * bmi", logit.link = TRUE)

dp <- addCondition(defB, dp, newvar = "int_B")
dp <- addCondition(defD, dp, newvar = "int_D")

Miss <- genMiss(dp, defM, 
  idvars = c("timeID","id", "period", "base_age", "age", "sex", "int_B", "int_D",
             "theta_bmi", "theta_dbp"))

dobs <- genObs(dp, Miss,
  idvars = c("timeID","id", "period", "base_age", "age", "sex", "int_B", "int_D",
             "theta_bmi", "theta_dbp"))
data<-rbind(dobs,data)
}

data.sim<-cbind(data,rep(1:n.sim,each=n*6))
names(data.sim)
data.sim<- data.sim[,-c(1,4,7:10)]
names(data.sim)[names(data.sim) == "V2"] <- "dataset"
nrow(data.sim)
summary(data.sim)
write.table(data.sim,"simulation.keith.long.1000.csv",sep=",",col.names=T,row.names=F)

A<-NULL
for (i in 1:n.sim){
B<-reshape(data.sim[((i-1)*nrow(data.sim[which(data.sim$dataset==i),])+1):((i-1)*nrow(data.sim[which(data.sim$dataset==i),])+nrow(data.sim[which(data.sim$dataset==i),])),], idvar = "id", timevar = "period", direction = "wide")
A<-rbind(A,B)
}
names(A)
summary(A)
nrow(A)

names(A)[names(A) == "age.0"] <- "age0"
names(A)[names(A) == "age.1"] <- "age1"
names(A)[names(A) == "age.2"] <- "age2"
names(A)[names(A) == "age.3"] <- "age3"
names(A)[names(A) == "age.4"] <- "age4"
names(A)[names(A) == "age.5"] <- "age5"
names(A)[names(A) == "bmi.0"] <- "bmi0"
names(A)[names(A) == "bmi.1"] <- "bmi1"
names(A)[names(A) == "bmi.2"] <- "bmi2"
names(A)[names(A) == "bmi.3"] <- "bmi3"
names(A)[names(A) == "bmi.4"] <- "bmi4"
names(A)[names(A) == "bmi.5"] <- "bmi5"
names(A)[names(A) == "dbp.0"] <- "dbp0"
names(A)[names(A) == "dbp.1"] <- "dbp1"
names(A)[names(A) == "dbp.2"] <- "dbp2"
names(A)[names(A) == "dbp.3"] <- "dbp3"
names(A)[names(A) == "dbp.4"] <- "dbp4"
names(A)[names(A) == "dbp.5"] <- "dbp5"
names(A)[names(A) == "sex.0"] <- "sex"
names(A)[names(A) == "dataset.0"] <- "dataset"
A<-A[,-c(8,11,13,16,18,21,23,26,28,31)]
nrow(A)
names(A)
write.table(A,"simulation.keith.wide.1000.csv",sep=",",col.names=T,row.names=F)

data1<-A[which(is.na(A$bmi0)==FALSE | is.na(A$bmi1)==FALSE | is.na(A$bmi2)==FALSE | is.na(A$bmi3)==FALSE | is.na(A$bmi4)==FALSE | is.na(A$bmi5)==FALSE),]
nrow(data1)

data2<-data1[which(is.na(data1$dbp0)==FALSE | is.na(data1$dbp1)==FALSE | is.na(data1$dbp2)==FALSE | is.na(data1$dbp3)==FALSE | is.na(data1$dbp4)==FALSE | is.na(data1$dbp5)==FALSE ),]
nrow(data2)

write.table(data2,"simulation.keith.1000.csv",sep=",",col.names=T,row.names=F)

end.simulation<-Sys.time()
time.simulation<- end.simulation-start.simulation
time.simulation

