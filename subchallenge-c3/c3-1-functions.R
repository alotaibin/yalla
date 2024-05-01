require(mvnfast)
require(this.path)
setwd(this.path::here())

# CDF, quantile function, Jacobian for transformation to delta-Laplace margins.

#------------------------------------------------------------------------------------------------------------------
# pdlaplace: univariate distribution function

pdlaplace<-function(z,mu,sigma,delta,lower.T=T)
{
  
  k=sqrt(gamma(1/delta)/gamma(3/delta))
  sigma=k*sigma
  nz<-length(z)
  result<-numeric(nz)
  if(length(delta)==1){delta<-rep(delta,nz)}
  if(lower.T==T){
    result[z<mu]<- 0.5*pgamma((((mu-z)/sigma)[z<mu])^delta[z<mu],shape=1/delta[z<mu],scale=1,lower.tail = F)
    result[z>=mu]<- 0.5+0.5*pgamma((((z-mu)/sigma)[z>=mu])^delta[z>=mu],shape=1/delta[z>=mu],scale=1)
  }else{
    
    result[z<mu]<- 1-0.5*pgamma((((mu-z)/sigma)[z<mu])^delta[z<mu],shape=1/delta[z<mu],scale=1,lower.tail = F)
    result[z>=mu]<- 0.5*pgamma((((z-mu)/sigma)[z>=mu])^delta[z>=mu],shape=1/delta[z>=mu],scale=1,lower.tail=F)
  }
  
  return(result)
}

#------------------------------------------------------------------------------------------------------------------
# qdlaplace: univariate quantile function

qdlaplace<-function(q,mu,sigma,delta)
{
  k=sqrt(gamma(1/delta)/gamma(3/delta))
  sigma=k*sigma
  nq<-length(q)
  result<-numeric(nq)
  if(length(mu)==1){mu<-rep(mu,nq)}
  if(length(sigma)==1){sigma<-rep(sigma,nq)}
  if(length(delta)==1){delta<-rep(delta,nq)}
  result[q<0.5]<- mu[q<0.5] - sigma[q<0.5]*qgamma(1-2*q[q<0.5],shape=1/delta[q<0.5],scale=1)^(1/delta[q<0.5])
  result[q>=0.5]<- mu[q>=0.5] + sigma[q>=0.5]*qgamma(2*q[q>=0.5]-1,shape=1/delta[q>=0.5],scale=1)^(1/delta[q>=0.5])
  return(result)
}


#------------------------------------------------------------------------------------------------------------------
# rdlaplace: univariate random number generation

rdlaplace<-function(n,mu,sigma,delta)
{
  k=sqrt(gamma(1/delta)/gamma(3/delta))
  sigma=k*sigma
  return(qdlaplace(runif(n),mu=mu,delta=delta,sigma=sigma))
}


#------------------------------------------------------------------------------------------------------------------
# ddlaplace: univariate density function

ddlaplace<-function(z,mu,sigma,delta,log=FALSE)
{
  k=sqrt(gamma(1/delta)/gamma(3/delta))
  sigma=k*sigma
  ld<- -abs((z-mu)/(sigma))^delta +log(delta)-log(2*sigma)-lgamma(1/delta)
  if(log==TRUE){return(ld)}
  else{return(exp(ld))}
}

#------------------------------------------------------------------------------------------------------------------
# dmvdlaplace: multivariate density function with Gaussian copula


# with faster MVN computation: requires matrix Sigma AND its cholesky factorization, chol(Sigma), which results in faster
# computation when used in apply()

dmvdlaplace<-function(z,mu,sigmad,SigmaChol,Sigma,delta,log=FALSE)
{
  zn<-qnorm(pdlaplace(z,mu=mu,sigma=sigmad,delta=delta,lower.T=F),mean=mu,sd=sqrt(diag(Sigma)),lower.tail =F)
  ld1<-dmvn(zn,mu=mu,sigma = SigmaChol,log=TRUE,isChol=TRUE)
  ld2<-sum(ddlaplace(z,mu=mu,sigma=sigmad,delta=delta,log=TRUE))-sum(dnorm(zn, mean=mu, sd=sqrt(diag(Sigma)), log=TRUE))
  if(log==TRUE){return(ld1+ld2)}
  else{return(exp(ld1+ld2))}
}

#------------------------------------------------------------------------------------------------------------------
# rmvdlaplace: multivariate random number generation with Gaussian copula

rmvdlaplace<-function(n,dim,mu,sigmad,Sigma,delta)
{
  if(dim(Sigma)[1]!=dim){stop("dim and Sigma do not match")}
  NU<-t(apply(rmvn(n,mu=rep(0,length=dim),sigma=Sigma),1,pnorm,sd=sqrt(diag(Sigma))))
  return(t(apply(NU,1,qdlaplace,mu=mu,sigma=sigmad,delta=delta)))
}

# Heffernan and Tawn nll

nll=function(par,X0,X1X2,covariates){
  
  alpha.mat=matrix(par[1:6], ncol=2, byrow=F )
  Y=as.matrix(cbind(rep(1,length(covariates[,1])),covariates)) 
  alpha.Y=tanh(Y%*%alpha.mat)  #t(Y%*%alpha.mat)
  beta.vec=par[7:8]
  sig.vec=par[9:10]
  mu.vec=par[11:12]
  delta.vec=par[13:14]
  rho=par[15]
  
  #Add in any constraints here
  if(any(alpha.Y< -1) | any(alpha.Y > 1 ) | any(beta.vec> 1 ) | any(beta.vec < 0) | any(sig.vec<0) | any(delta.vec < 0) )  return(1e10)
  
  cov.mat=diag(sig.vec)%*%matrix(c(1,rho,rho,1),2,2)%*%diag(sig.vec)
  
  nll<-apply(cbind(X0,X1X2,alpha.Y),1,function(x){
    mean=x[4:5]*x[1]+(x[1]^beta.vec)*mu.vec
    var=diag(c(x[1]^beta.vec))%*%cov.mat%*%diag(c(x[1]^beta.vec))
    dmvdlaplace(x[2:3]-mean,mu= c(0,0),sigmad=sqrt(diag(var)),SigmaChol=chol(var),Sigma=var,
                delta=delta.vec,log=T)
  })
  
  
  if(is.finite(sum(nll))){
    return(-sum(nll))
  }else{return(1e10)}
  
  
}

##### Upload data:
Coputopia <- as.data.frame(read_csv("../Data/Coputopia.csv",show_col_types = FALSE))
Coputopia$Season = as.numeric(as.factor(Coputopia$Season)) #(S1,S2) as (1,0)
data=as.matrix(Coputopia[,3:5])
Y1 <- data[,1]
Y2 <- data[,2]
Y3 <- data[,3]
##### Transform Y~Gumbel(0,1) to U~uniform(0,1):
U1 <- exp(-exp(-X1))
U2 <- exp(-exp(-X2))
U3 <- exp(-exp(-X3))
##### Transform X's to Laplace distribution:
LaplaceTrans <- function(u,location,scale){
  return(location - sign(u - 0.5) * scale * (log(2) + ifelse(u < 0.5, log(u), log1p(-u))))
}
L1 <- LaplaceTrans(U1,0,1)
L2 <- LaplaceTrans(U2,0,1)
L3 <- LaplaceTrans(U3,0,1)
Y=as.matrix(cbind(L1,L2,L3))
covariates=as.matrix(Coputopia[,1:2])

###########################################  
##########  Model 1: Y2Y3 | Y1   ##########  
Y1=Y[,1]
Y2Y3=Y[,c(2,3)]
u=LaplaceTrans(0.95,0,1)
# Identify the values that exceed the threshold 'u'
Y1.exceed=Y1[Y1>u]
Y2Y3.exceed=Y2Y3[Y1>u,]

# Plot the bivariate data
par(mfrow=c(1,1))
plot(Y1,Y2Y3[,1],pch=20,cex=.9)
points(Y1.exceed,Y2Y3.exceed[,1],pch=20,cex=.9,col="red")

# set initial parameter values for the optimization function 
init.par=rep(0.5,15)
covariates1=as.matrix(cbind(Coputopia[Y1>u,1],(Coputopia[Y1>u,2])^3))

# fitting the HT model:
start_time <- Sys.time()
fit1<-optim(par=init.par,fn=nll,X0=Y1.exceed,X1X2=Y2Y3.exceed,covariates=covariates1) 
end_time <- Sys.time()
Tot.Time <- end_time - start_time

# parameter estimates & AIC value for model 1
fit1$par
AIC1 = 2*(15+2) - 2*fit1$value

###########################################
##########  Model 2: Y1Y3 | Y2   ##########  
Y2=Y[,2]
Y1Y3=Y[,c(1,3)]
# Identify the values that exceed the threshold 'u'
Y2.exceed=Y2[Y2>u]
Y1Y3.exceed=Y1Y3[Y2>u,]

# Plot the bivariate data
par(mfrow=c(1,1))
plot(Y2,Y1Y3[,1],pch=20,cex=.9)
points(Y2.exceed,Y1Y3.exceed[,1],pch=20,cex=.9,col="red")

# set initial parameter values for the optimization function
init.par=rep(0.5,15)
covariates2=as.matrix(cbind(Coputopia[Y2>u,1],(Coputopia[Y2>u,2])^3))

# fitting the HT model:
start_time <- Sys.time()
fit2<-optim(par=init.par,fn=nll,X0=Y2.exceed,X1X2=Y1Y3.exceed,covariates=covariates2) 
end_time <- Sys.time()
Tot.Time <- end_time - start_time

# parameter estimates & AIC value for model 2
fit2$par
AIC2 = 2*(15+2) - 2*fit2$value

###########################################  
##########  Model 3: Y1Y2 | Y3   ##########  
Y3=Y[,3]
Y1Y2=Y[,1:2]
# Identify the values that exceed the threshold 'u'
Y3.exceed=Y3[Y3>u]
Y1Y2.exceed=Y1Y2[Y3>u,]

# Plot the bivariate data
par(mfrow=c(1,1))
plot(Y3,Y1Y2[,1],pch=20,cex=.9)
points(Y3.exceed,Y1Y2.exceed[,1],pch=20,cex=.9,col="red")

# set initial parameter values for the optimization function. 
init.par=rep(0.5,15)
covariates3=as.matrix(cbind(Coputopia[Y3>u,1],(Coputopia[Y3>u,2])^3))

# fitting the HT model:
start_time <- Sys.time()
fit3<-optim(par=init.par,fn=nll,X0=Y3.exceed,X1X2=Y1Y2.exceed,covariates=covariates3) 
end_time <- Sys.time()
Tot.Time <- end_time - start_time

# parameter estimates & AIC value for model 3
fit3$par
AIC3 = 2*(15+2) - 2*fit3$value



##############################################################
#############  Parameter Estimates for all models  ###########
##############################################################

alpha.mat1=matrix(fit1$par[1:6], ncol=2, byrow=F )
Y11=as.matrix(cbind(rep(1,length(covariates1[,1])),covariates1)) 
alpha.Y1=tanh(Y11%*%alpha.mat1 )
beta.vec1 <- fit1$par[7:8]

alpha.mat2=matrix(fit2$par[1:6], ncol=2, byrow=F )
Y22=as.matrix(cbind(rep(1,length(covariates2[,1])),covariates2)) 
alpha.Y2=tanh(Y22%*%alpha.mat2 ) 
beta.vec2 <- fit2$par[7:8]

alpha.mat3=matrix(fit3$par[1:6], ncol=2, byrow=F )
Y33=as.matrix(cbind(rep(1,length(covariates3[,1])),covariates3)) 
alpha.Y3=tanh(Y33%*%alpha.mat3)  
beta.vec3 <- fit3$par[7:8]

#################################################
############# Residuals for all models ##########
#################################################

##residuals = {(Y2Y3) - (a1a2)Y1} / Y1^(b1b2)
residuals1 = (Y1Y2.exceed - alpha.Y1*Y1.exceed)/(cbind(Y1.exceed^beta.vec1[1],Y1.exceed^beta.vec1[2]))  
residuals2 = (Y1Y2.exceed - alpha.Y2*Y2.exceed)/(cbind(Y2.exceed^beta.vec2[1],Y2.exceed^beta.vec2[2]))  
residuals3 = (Y1Y2.exceed - alpha.Y3*Y3.exceed)/(cbind(Y3.exceed^beta.vec3[1],Y3.exceed^beta.vec3[2]))  


Model1 <- data.frame(Model1_a1=alpha.Y1[,1], Model1_a2=alpha.Y1[,2], res_z1=residuals1[,1],res_z2=residuals1[,2])
Model2 <- data.frame(Model2_a1=alpha.Y2[,1], Model2_a2=alpha.Y2[,2], res_z1=residuals2[,1],res_z2=residuals2[,2])
Model3 <- data.frame(Model3_a1=alpha.Y3[,1], Model3_a2=alpha.Y3[,2], res_z1=residuals3[,1],res_z2=residuals3[,2])

beta <- data.frame(Model1 = beta.vec1, Model2 = beta.vec2, Model3 = beta.vec3, row.names = c("b1","b2"))

write.csv(Model1,'Model1_alpha_res.csv', row.names = FALSE)
write.csv(Model2,'Model2_alpha_res.csv', row.names = FALSE)
write.csv(Model3,'Model3_alpha_res.csv', row.names = FALSE)
write.csv(beta,'all_beta.csv')
