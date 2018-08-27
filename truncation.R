# Truncated regression

library(magrittr);
library(dplyr)

# Simulated data

set.seed(1056)                      # set seed to replicate example
nobs = 2500                         # number of obs in model
x1 <- runif(nobs,-5,5)               # random uniform variable
alpha = 10                        # intercept
beta = 7                            # angular coefficient
beta2 = 0.5
beta3 = -0.75
xb <- alpha + beta*x1 + beta2*x1^2  + beta3*x1^3            # linear predictor, xb
sd <- 1.5                             # Standard deviation
y <- rnorm(nobs, xb, sd = sd)       # create y as  random normal variate
cc <- function(x){5*(1+0.3*log(x^2))}
regdat <- data.frame(x1,y)
plot(regdat)


# Truncated y
truncdat <- regdat %>% filter(.,y >=cc(x1))
ntrunc <- nrow(truncdat)
ytrunc <- truncdat$y
xtrunc <- truncdat$x1
plot(truncdat)

# Likelihood
trunc_likelihood <- function(par){
  alpha = par[1]
  beta = par[2]
  beta2 = par[3]
  beta3 = par[4]
  sigma = par[5]
  
  xb =  alpha + beta*xtrunc + beta2*xtrunc^2 + beta3*xtrunc^3 

  lnL <-  -log(sigma) + dnorm((ytrunc - xb)/sigma, log = TRUE) - pnorm((xb - cc(xtrunc))/sigma, log.p = TRUE)
  LL <- sum(lnL)
  return(LL)
}

# Prior
low <- c(rep(-50,4),1e-5)
up <- c(rep(50,4),10)

prior <- createUniformPrior(lower = low,
                            upper = up)


setup <- createBayesianSetup(likelihood = trunc_likelihood,prior = prior)

settings <- list(iterations = 1e5,adaptation = 0.25,
                 burnin = 2e4, message = T,nrChains = 1)


system.time(
  res <- runMCMC(bayesianSetup = setup, settings = settings,sampler = "DREAMzs")
)

summary(res)


codaObject = getSample(res, start = 1E3, coda = TRUE)

getmcmc_var <- function(outjags=outjags, vars = vars){
  as.data.frame(do.call(rbind, outjags[,vars]))
}
ss <- getmcmc_var(codaObject,vars = c("par 1","par 2","par 3","par 4","par 5"))
colnames(ss) <- c("alpha","beta","beta2","beta3","sd")
index <- sample(seq(1:nrow(ss)),250,replace=FALSE)
ss <- ss[index,]

xpred <- seq(min(x1),max(x1),length.out = 250)
df <- NULL
for(i in 1:250){
  temp_df <- data.frame(x=xpred,y= (ss$alpha[i] + ss$beta[1]*xpred + ss$beta2[1]*xpred^2 +
  ss$beta3[1]*xpred^3),col=rep(i:i, each=250))
df <- rbind(df,temp_df)  
}


ggplot(data=truncdat,aes(x=x1,y=y)) +
  
  geom_point() +
  #  geom_segment(data = filter(censdat,y==L_limit),
  #               mapping=aes(x=x1, y=y, xend=x1, yend = y-5),size=0.1,
  #               colour=cens,arrow = arrow()) +
  geom_point(data=filter(regdat,y <= cc(x1)),mapping=aes(x=x1,y=y),color="gray50")+
  coord_cartesian(ylim=c(-5,50)) +
  theme_pubr() +
  scale_color_stata() +
  scale_shape_stata()+
  xlab("x") + ylab("y") +
  stat_smooth(formula=y ~ poly(x1, 3, raw=TRUE),linetype="dashed",colour="red",se=F) +
  geom_line(data=df,aes(x = x, y = y,group=col),
              alpha = 0.1, color = "green",size=0.3) 






#+
#  geom_abline(slope = mean(ss$beta), 
#              intercept = mean(ss$alpha), 
#              color = "orange3", size = 1) 






