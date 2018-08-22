# Truncated regression

library(magrittr);
library(dplyr)

# Simulated data

set.seed(1056)                      # set seed to replicate example
nobs = 1000                         # number of obs in model
x1 <- runif(nobs,0,5)               # random uniform variable
alpha = 4                         # intercept
beta = 7                            # angular coefficient
xb <- alpha + beta*x1               # linear predictor, xb
sd <- 1.5                             # Standard deviation
y <- rnorm(nobs, xb, sd = sd)       # create y as  random normal variate
c <- 20
regdat <- data.frame(x1,y)
# Truncated y


truncdat <- regdat %>% filter(.,y >=c)

ytrunc <- truncdat$y
xtrunc <- truncdat$x1


trunc_likelihood <- function(par){
  alpha = par[1]
  beta = par[2]
  sigma = par[3]
  xb =  alpha + beta*xtrunc
  lly  =  sum(log(pnorm((c - xb)/sigma)))
  llobs = sum(dnorm((ytrunc - xb)/sigma,sd = sigma,log = T))
  return(llobs + lly)
}


low <- c(-50,-50,1e-5)
up <- c(50,50,100)

prior <- createUniformPrior(lower = low,
                            upper = up)
setup <- createBayesianSetup(likelihood = trunc_likelihood,prior = prior)

settings <- list(iterations = 5e4,adaptation = 0.25,
                 burnin = 1e4, message = T,nrChains = 1)


system.time(
  res <- runMCMC(bayesianSetup = setup, settings = settings,sampler = "DREAMzs")
)

summary(res)


codaObject = getSample(res, start = 500, coda = TRUE)

getmcmc_var <- function(outjags=outjags, vars = vars){
  as.data.frame(do.call(rbind, outjags[,vars]))
}
ss <- getmcmc_var(codaObject,vars = c("par 1","par 2","par 3"))
colnames(ss) <- c("alpha","beta","sd")
index <- sample(seq(1:nrow(ss)),100,replace=FALSE)
ss <- ss[index,]




ggplot(data=truncdat,aes(x=x1,y=y)) +
  
  geom_point() +
  #  geom_segment(data = filter(censdat,y==L_limit),
  #               mapping=aes(x=x1, y=y, xend=x1, yend = y-5),size=0.1,
  #               colour=cens,arrow = arrow()) +
  geom_point(data=filter(regdat,y <= c),mapping=aes(x=x1,y=y),color="gray90")+
  coord_cartesian(ylim=c(-2,50)) +
  theme_pubr() +
  scale_color_stata() +
  scale_shape_stata()+
  xlab("x") + ylab("y") +
  geom_smooth(method="lm",linetype="dashed",se=F) +
  geom_abline(data=ss,aes(intercept = alpha, slope = beta),
              alpha = 0.1, color = "green",size=0.3) +
  geom_abline(slope = mean(ss$beta), 
              intercept = mean(ss$alpha), 
              color = "orange3", size = 1) 
