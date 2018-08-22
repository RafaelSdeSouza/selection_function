# Censoring regression

library(magrittr);
library(dplyr)

# Simulated data

set.seed(1056)                      # set seed to replicate example
nobs = 500                         # number of obs in model
x1 <- runif(nobs,0,5)               # random uniform variable
alpha = 2.5                         # intercept
beta = 7                            # angular coefficient
xb <- alpha + beta*x1               # linear predictor, xb
sd <- 2                             # Standard deviation
y <- rnorm(nobs, xb, sd = sd)       # create y as  random normal variate
L_limit <- 20
regdat <- data.frame(x1,y)
# Censored y


censdat <- regdat %>% mutate(.,cens = ifelse(y <= L_limit, "yes","no")) %>%
mutate(y = replace(y,y <= L_limit, L_limit))

ycens <- censdat$y
xcens <- censdat$x1


tobit_likelihood <- function(par){
  alpha = par[1]
  beta = par[2]
  sigma = par[3]
  Id = ycens > L_limit  
  mu  =  alpha + beta*xcens
  lly  = sum((1 - Id)*pnorm((L_limit - mu),sd=sigma, log.p = T))
  llobs = sum(Id*dnorm(ycens,mean = mu ,sd = sigma,log = T))
  return(llobs + lly)
}


low <- c(-50,-50,1e-5)
up <- c(50,50,100)

prior <- createUniformPrior(lower = low,
                            upper = up)
setup <- createBayesianSetup(likelihood = tobit_likelihood,prior = prior)

settings <- list(iterations = 5e4,adaptation = 0.25,
                 burnin = 1e4, message = T,nrChains = 1)


system.time(
  res <- runMCMC(bayesianSetup = setup, settings = settings,sampler = "DREAMzs")
)

codaObject = getSample(res, start = 500, coda = TRUE)

getmcmc_var <- function(outjags=outjags, vars = vars){
  as.data.frame(do.call(rbind, outjags[,vars]))
}
ss <- getmcmc_var(codaObject,vars = c("par 1","par 2","par 3"))
colnames(ss) <- c("alpha","beta","sd")
index <- sample(seq(1:nrow(ss)),100,replace=FALSE)
ss <- ss[index,]


summary(mod <- lm(ycens ~ xcens))    



ggplot(data=censdat,aes(x=x1,y=y)) +
  
  geom_point(aes(group=cens,shape=cens,colour=cens)) +
#  geom_segment(data = filter(censdat,y==L_limit),
#               mapping=aes(x=x1, y=y, xend=x1, yend = y-5),size=0.1,
#               colour=cens,arrow = arrow()) +
  geom_point(data=filter(regdat,y <= L_limit),mapping=aes(x=x1,y=y),color="gray90")+
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
  











# Fit with R function lm
summary(mod <- lm(ycens ~ xcens))          # model of the synthetic data.

# Plot Output
ypred <- predict(mod,type="response")                # prediction from the model
plot(xcens,ycens,pch = 19,col="red")                          # plot scatter
lines(xcens,ypred,col = 'grey40',lwd=2)                   # plot regression line
segments(xcens,fitted(mod),xcens,ycens,lwd = 1,col="gray70")     # add the residuals





# Fit uncensored data with R function lm
summary(mod2 <- lm(y ~ x1,data = censdat))          # model of the synthetic data.






library(VGAM) # load package


fit1 <- vglm(censdat$y ~ censdat$x1,tobit(Lower = 10))
coef(fit1, matrix = TRUE)
summary(fit1)

# Plot Output
ypred_cens <- predict(mod2,type="response")                # prediction from the model
plot(censdat$x1,censdat$y,pch=19,col="red")                          # plot scatter
lines(x1,ypred_cens,col='grey40',lwd=2,data=censdat)                   # plot regression line
segments(x1,fitted(mod2),x1,y,lwd=1,col="gray70",data=censdat)     # add the residuals
