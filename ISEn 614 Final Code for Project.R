library(gdata)
library(readxl)
library(RSpectra)
library(EnvStats)
library(matlib)
getwd()
p=read_xlsx("Project_dataset.xlsx", col_names= FALSE)
data_set=p
summary(data_set)
dim(data_set)
#data_set is the original data we have
#checking missing values
miss_vals=c()
for(x in 1:ncol(data_set)){
  p=sum(is.na(data_set[,x]))
  miss_vals=c(miss_vals,p)
}
miss_vals
#Percentage of variance explained:
var_perc=sum(res$values[1:3])/sum(res$values)
var_perc
#---------------------
#This set of code just shows that all the data attributes follows the normal distribution:
cor(data_set) #Calculates the correlations, can be helpful in report
par(mfrow=c(2,5))
x=data_set1
for(i in 1:ncol(data_set1)){
  m=mean(data_set1[,i])
  s=sd(data_set1[,i])
  hist(data_set1[,i], freq=FALSE, breaks=8, col="White")
  curve(dnorm(x, mean=m, sd=s),col="Red",add=TRUE)
}
m=mean(data_set1[,1])
s=sd(data_set1[,1])
hist(data_set1[,1], freq=FALSE, breaks=8, col="White")
curve(dnorm(x, mean=m, sd=s),col="Red",add=TRUE)
#--------------------
#--------------------------------
#The PCA Part:
#PCA:
p_cov=cov(data_set)
dim(p_cov)
res=eigs(p_cov,209,which="LM")
eig_vals=res$values
U=res$vectors
#plotting the eigen values"
K_vals=res$values[c(1:4)]
cheg=c()
perc_var=sum(res$values[c(1:4)])/sum(res$values)
perc_var
#Generating the data set from PCA:
U1=U[,c(1:4)]
k=c()
for(i in 1:nrow(data_set)){
  x=as.matrix(data_set[i,])
  z=x%*%U1
  k=as.data.frame(rbind(k,z))
}
dim(k)
#Data set for T^2 chart:
#Covariance and its inverse
t_cov=cov(k)
dim(t_cov)
inv_cov=inv(t_cov)
#mean vector computation
l1=c()
for(i in 1:ncol(k)){
  j=k[,i]
  l=mean(j)
  l1=c(l1,l)
}
o1=k[1,]
j=as.matrix(o1-l1)
dim(j)
#----------------------------------
t=c()
for(i in 1:nrow(k)){
  o1=k[i,]
  p1=as.matrix(o1-l1)
  t_sq=p1%*%inv_cov%*%t(p1)
  t=c(t,t_sq)
}
t
length(t)
#So t is the 552 t^2 statistics, now we just need to plot them in chart
#-------------------------------------------
#potting
par(mfrow=c(1,1))
plot(t)
abline(a=9.48,b=0) #a=UCL here:
#-------------------------------------------------
#Main Function for t^2 statistic and iterative process of removing out of control data points:

t2_stat=function(t,lu,ucl){
  r_pt=which(t>ucl)
  #Covariance and its inverse
  k1=lu[-r_pt,]
  t_cov=cov(k1)
  inv_cov=inv(t_cov)
  #mean vector computation
  l1=c()
  for(i in 1:ncol(k1)){
    j=k1[,i]
    l=mean(j)
    l1=c(l1,l)
  }
  t1=c()
  for(i in 1:nrow(k1)){
    o1=k1[i,]
    p1=as.matrix(o1-l1)
    t_sq=p1%*%inv_cov%*%t(p1)
    t1=c(t1,t_sq)
  }
  plot(t1)
  abline(a=ucl,b=0)
  return(list(t1,k1))
}
length(t)
dim(k)
i=0
j=1
#This one loop will do all the analysis for you, i just automated the whole process:
for(i in 1:100){
  if(any(t>9.48)==TRUE){
    newt=t2_stat(t,k,9.48)
    i=0
    t=newt[[1]]
    k=newt[[2]]
    j=j+1
  }else{
    t1=newt[[1]]
    k1=newt[[2]]
    break
  }
  
}
j #gives you the number of iterations it took to remove all the points.
t1 #The final t^2 points which are all in control
k1 #The final data_set containing 4 PCS which is fully in control
#-----------------------------
#--------------------------------
#MDL Code:
eig_vals
eig_vals1=sort(eig_vals)
eig_vals1
h=length(eig_vals1)-1
h
kr=c()
for(i in 0:h){
  e1=eig_vals1[209-i]
  ar=mean(e1)
  gr=geoMean(e1)
  MDL=(552*(209-i)*log(ar/gr))+(i*(418-i)*log(552)/2)
  kr=c(kr,MDL)
}
kr
which.min(kr)
#------------------------------
#Scree Plot for first 10 Eigen values:
lines(eig_vals[1:10])
#-----------------------------------------------------------