
# install.packages('survival')
# install.packages('NMOF')
# install.packages('ncvreg')
# install.packages('SIS')

library(survival)
library(NMOF)
library(ncvreg)
library(SIS)

### plot K-M plot
library(survival)
# install.packages("ggsurvplot")
library(ggsurvfit)
# install.packages("gtsummary")
# install.packages("survminer")
library('survminer')
library(gtsummary)

getwd()
setwd("C:/Users/xxx/xxx/FYP2_code")

system("R CMD SHLIB strcox-dh.c")
dyn.load('TSMCPCOX.dll')

##########################################################################################
### before running the following code, please load all functions in the function file 
##########################################################################################

peacedata = read.csv("peacedata.csv",header=TRUE)
head(peacedata)
total = dim(peacedata)[1]

x.trt = peacedata$tx
x.age = peacedata$age
x.bp = peacedata$sysbp
x.gender = peacedata$gender
x.diab = peacedata$hidiabet
x.hyper = peacedata$hihypert

X = cbind(x.trt, x.age, x.bp, x.gender, x.diab,x.hyper)
X = as.data.frame(X)
Z = peacedata$age
Y = peacedata$t2death
status = peacedata$death

# conduct subgroup identification
output= tsmcpcox(Y = Y, Delta=status ,X1=X,Z=Z, method = "MCP")
output


# save(output,file = "changepoint_save.Rdata")
change_pt = output$change.points
group1 = peacedata[peacedata$age<=67,]
d1=dim(group1)[1]

group2 = peacedata[peacedata$age>67 & peacedata$age<=70,]
d2=dim(group2)[1]

group3 = peacedata[peacedata$age>70,]
d3=dim(group3)[1]

cox_group1 = coxph(Surv(t2death,death)~., ties = 'breslow', data = group1)
summary(cox_group1)

cox_group2 = coxph(Surv(t2death,death)~., ties = 'breslow', data = group2)
summary(cox_group2)

cox_group3 = coxph(Surv(t2death,death)~., ties = 'breslow', data = group3)
summary(cox_group3)

############################################
## model comparison (with classical coxph)##
############################################

coxall = coxph(Surv(t2death,death)~., ties = 'breslow', data = peacedata)
summary(coxall)

# coxall model selection
cox_inhomo_1<-coxph(Surv(t2death,death)~tx*age*sysbp*gender*hidiabet*hihypert, ties="exact",data = peacedata)
stepwise_cox<-step(coxall,list(lower=~1,upper=formula(cox_inhomo_1)), direction = "both")
summary(stepwise_cox)
cox_step = coxph(formula = Surv(t2death, death) ~ tx + age + gender + hidiabet + 
                   hihypert + age:hidiabet + age:hihypert + tx:gender, 
                 data = peacedata, ties = "breslow")
summary(cox_step)
BIC(cox_step)

## model comparision withr respect to BIC
# BIC.submodels = BIC(cox_group1)+BIC(cox_group2)+BIC(cox_group3)
BIC.full = BIC(coxall)

l1 = cox_group1$loglik[2]
l2 = cox_group2$loglik[2]
l3 = cox_group3$loglik[2]

class(l1)
n = d1+d2+d3
k = 6*3
BIC.submodels = k*log(n)-2*(l1+l2+l3)

BIC.full
BIC(cox_step)


#####################################
###  KM- plot##
#####################################

data1 = peacedata
data1$agegroup = 0
data1$agegroup[data1$age<=67] = "group1"
data1$agegroup[data1$age>67&data1$age<=70] = "group2"
data1$agegroup[peacedata$age>70] = "group3"
data1$agegroup = as.factor(data1$agegroup)
# data1$agegroup
head(data1)

fit1 = survfit(Surv(t2death, death) ~ agegroup, data = data1)
ggsurvplot(fit1, conf.int = TRUE,pval=TRUE, ggtheme = theme_bw(), data = data1)

