# Load required libraries
require(gridExtra)
require(ggplot2)
require(rmutil)
require(readr)
require(this.path)
setwd(this.path::here())

# Calculate the 95th percentile of the Laplace distribution
u=qlaplace(0.95)

# Read model parameters and residuals from CSV files
betas <- read.csv("../Data/all_beta.csv")
input1 <- read.csv("../Data/Model1_alpha_res.csv")
input2 <- read.csv("../Data/Model2_alpha_res.csv")
input3 <- read.csv("../Data/Model3_alpha_res.csv")

# Store model parameters in a list structure for easy access
pars=vector("list",3)
pars[[1]]=list("beta"=t(as.matrix(betas$Model1)),"alpha.vec"=cbind(input1$Model1_a1,input1$Model1_a2),"resid"=cbind(input1$res_z1,input1$res_z2))
pars[[2]]=list("beta"=t(as.matrix(betas$Model2)),"alpha.vec"=cbind(input2$Model2_a1,input2$Model2_a2),"resid"=cbind(input2$res_z1,input2$res_z2))
pars[[3]]=list("beta"=t(as.matrix(betas$Model3)),"alpha.vec"=cbind(input3$Model3_a1,input3$Model2_a2),"resid"=cbind(input3$res_z1,input3$res_z2))

# Define a function to simulate data based on input parameters and conditions
rsim=function(beta,alphas,resid,v,cond.ind){
  
  X0=v+rexp(1)
  
  alpha.ind=sample(1:nrow(alphas),1)
  resid.ind=sample(1:nrow(resid),1)
  
  X1X2=(alphas[alpha.ind,])*X0+X0^beta*resid[resid.ind,]
  
  if(cond.ind==1) return(c(X0,X1X2))
  if(cond.ind==2) return(c(X1X2[1],X0,X1X2[2]))
  if(cond.ind==3) return(c(X1X2,X0))
  
}

rsim_all=function(pars,n,B,v){
  
  cond.inds=as.matrix(sample(1:3,B,replace=T))
  
  out<-t(apply(cond.inds,1,function(x){
    rsim(pars[[x]]$beta,pars[[x]]$alpha.vec,pars[[x]]$resid,v,x)
    
  }))
  
  import.prob=apply(out,1,function(x){
    1/sum(x>v)
  })
  import.prob = import.prob/sum(import.prob)
  
  samp.inds = sample(1:length(import.prob),size=n,prob=import.prob)
  
  return(out[samp.inds,])
}

# Read the data for C3
data <- read.csv("../Data/Coputopia.csv")
X=as.matrix(data[,3:5])

# Transform the data
X_L = qlaplace(exp(-exp(-X))) 
maxs = apply(X_L, 1, max)  # find the max value per row
max.exceed.inds = which(maxs > u)  # indices of max exceeding threshold
max.nonexceed.inds = which(maxs <= u)  # indices of max not exceeding threshold


n=10^7 # Number of simulations

start.time <- Sys.time()
boo<-rbinom(n=n,size=1,prob=length(max.exceed.inds)/nrow(X_L)  ) 
sims.exceeds<-rsim_all(pars,n=sum(boo),B=5*sum(boo),v=u)
end.time <- Sys.time() 
Tot.time <- end.time - start.time

sims=matrix(nrow=n,ncol=3)
sims[1:sum(boo),]=sims.exceeds
orig.inds=sample(max.nonexceed.inds,size=sum(boo==0),replace=T)
sims[-(1:sum(boo)),]=X_L[orig.inds,]

# Back-transform simulated data
sims=-log(-log(plaplace(sims)))


## save the simulations
write.csv(sims,"sims.results.csv",row.names = F)

###############################################################################
## Compute the probabilities for C3:
y=6
v=7
m=-log(log(2))

## Prob 1
P11 <- mean(sims[,1] > y & sims[,2] > y & sims[,3] > y) #0.0001933

## Prob 2
P22 <- mean(sims[,1] > v & sims[,2] > v & sims[,3] < m)  #1.71e-05
AnswerC3 <- c(P11, P22)

## Save the results for C3:
write.csv(AnswerC3,"AnswerC3.csv", row.names = F)


###########################################################
############# Diagnostics plots for Model: QQ plots  ######
###########################################################
#input.sims <- read.csv("sims.results.csv")
input.sims <- sims
min.sims <- apply(input.sims, 1, min)
max.sims <- apply(input.sims, 1, max)
sum.sims <- apply(input.sims, 1, sum)

min.data <- apply(data[3:5], 1, min)
max.data <- apply(data[3:5], 1, max)
sum.data <- apply(data[3:5], 1, sum)

### QQ plots
par(mfrow=c(1,3), mar = c(4,4,4,4))
qqplot(min.sims, min.data, xlab = "Simulated with u=0.95", ylab = "Observed min", ylim=c(-4,15), xlim=c(-4,15), main = "", asp=1)
abline(0, 1, col = "lightgrey")

qqplot(max.sims, max.data, xlab = "Simulated with u=0.95", ylab = "Observed max", ylim=c(-2,25), xlim=c(-2,25), main = "", asp=1)
abline(0, 1, col = "lightgrey")

qqplot(sum.sims, sum.data, xlab = "Simulated with u=0.95", ylab = "Observed sum", ylim=c(-6,25), xlim=c(-6,25), main = "",asp=1)
abline(0, 1, col = "lightgrey")

##########################################################
# Calculate quantiles for x and y to be ploted with ggplot

min.sims1 <- quantile(min.sims, probs = seq(0, 1, length.out = length(min.sims)))
min.data1 <- quantile(min.data, probs = seq(0, 1, length.out = length(min.data)))
max.sims1 <- quantile(max.sims, probs = seq(0, 1, 0.01))
max.data1 <- quantile(max.data, probs = seq(0, 1, 0.01))
sum.sims1 <- quantile(sum.sims, probs = seq(0, 1, 0.01))
sum.data1 <- quantile(sum.data, probs = seq(0, 1, 0.01))


min_data <- data.frame(x=min.sims1,y=min.data1)
max_data <- data.frame(x=max.sims1,y=max.data1)
sum_data <- data.frame(x=sum.sims1,y=sum.data1)

q1 <- ggplot(min_data, aes(x = x, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.2) +
  coord_fixed(ratio = 1) +
  labs(title = "",
       x = "Observed min",
       y = "Theoretical min")

q2 <- ggplot(max_data, aes(x = x, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed" , color = "red", size = 0.2) +
  coord_fixed(ratio = 1) +
  labs(title = "",
       x = "Observed max",
       y = "Theoretical max")

q3 <- ggplot(sum_data, aes(x = x, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.2) +
  coord_fixed(ratio = 1) +
  labs(title = "",
       x = "Observed sum",
       y = "Theoretical sum")

combined_plot <- grid.arrange(q1, qq, q3, ncol = 3) 



