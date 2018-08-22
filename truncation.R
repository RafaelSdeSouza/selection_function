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
c <- 5
regdat <- data.frame(x1,y)
plot(regdat)
# Truncated y


truncdat <- regdat %>% filter(.,y >=c)
ntrunc <- nrow(truncdat)
ytrunc <- truncdat$y
xtrunc <- truncdat$x1


trunc_likelihood <- function(par){
  alpha = par[1]
  beta = par[2]
  beta2 = par[3]
  beta3 = par[4]
  sigma = par[5]
  
  xb =  alpha + beta*xtrunc + beta2*xtrunc^2 + beta3*xtrunc^3 

  lnL <-  -log(sigma) + dnorm((ytrunc - xb)/sigma, log = TRUE) - pnorm((xb - c)/sigma, log.p = TRUE)
  LL <- sum(lnL)
  return(LL)
}


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
  geom_point(data=filter(regdat,y <= c),mapping=aes(x=x1,y=y),color="gray50")+
  coord_cartesian(ylim=c(-2,50)) +
  theme_pubr() +
  scale_color_stata() +
  scale_shape_stata()+
  xlab("x") + ylab("y") +
  geom_smooth(method="loess",linetype="dashed",colour="red",se=F) +
  geom_line(data=df,aes(x = x, y = y,group=col),
              alpha = 0.1, color = "green",size=0.3) 
#+
#  geom_abline(slope = mean(ss$beta), 
#              intercept = mean(ss$alpha), 
#              color = "orange3", size = 1) 






# Construct data dictionary
X <- model.matrix(~ 1 + xtrunc)
K <- ncol(X)
model.data <- list(Y = ytrunc,          # Response variable
                   X = X,               # Predictors
                   K = K,               # Number of predictors including the intercept
                   N = ntrunc,             # Sample size
                   c = c
)

# Model set up
NORM <- "model{
# Diffuse normal priors for predictors
for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001) }
# Uniform prior for standard deviation
tau <- pow(sigma, -2) # precision
sigma ~ dunif(0, 100) # standard deviation
# Likelihood function
for (i in 1:N){
Y[i] ~ dnorm(mu[i],tau)T(c,)
mu[i] <- eta[i]
eta[i] <- inprod(beta[], X[i,])
}
}"

# Initial values
inits <- function () {
  list(beta = rnorm(K, 0, 0.01))
}

# Parameters to be displayed
params <- c("beta", "sigma")

# MCMC
normfit <- jags(data = model.data,
                inits = inits,
                parameters = params,
                model = textConnection(NORM),
                n.chains = 3,
                n.iter = 15000,
                n.thin = 1,
                n.burnin = 10000)

print(normfit, intervals = c(0.025, 0.975), digits = 2)


