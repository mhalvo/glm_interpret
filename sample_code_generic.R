##########################################
# sample_code_generic.R                  #
#    Code used to generate figures       #
#     and table estimates for simulated  #
#     paper example: meant to be adapted #
#     for individual use cases           #
##########################################
# Max Halvorson #
#   2020-05-07  #
#################

rm(list=ls())
library('MASS')
library('ggplot2')
options(max.print=2000, scipen=7, tibble.print_min=70)

setwd('C:/Users/Max Halvorson/Dropbox/GLM Main Effects Paper')

##########################
# Simulated data example #
##########################
### Simulate data
set.seed(1337) 
b0 <- 0        # Intercept
b1 <- -2       # X1 Effect
b2 <- 1        # X2 Effect
n <- 500       # Sample Size
mu <- rep(0,2) # Specify means
S <- matrix(c(  1,-.1,
              -.1,  1),nrow=2,ncol=2)   # Specify covariance matrix
sigma <- 1                              # Level 1 error
rawvars <- mvrnorm(n=n, mu=mu, Sigma=S) # Simulate predictors from a multivariate normal distribution

id <- seq(1:n)
eij <- rep(rnorm(id, 0, sigma))         # Add normally distributed error
xb <- (b0) + (b1)*rawvars[,1] + (b2)*rawvars[,2] + eij

pr <- 1 / (1 + exp(-xb))        # pass through logistic function to convert xb to probabilities
y <- rbinom(n, 1, pr)           # bernoulli response variable
df <- data.frame(id=id,
                 y=y,
                 x1=rawvars[,1],
                 x2=rawvars[,2])

##################################
# Analysis for graphs and tables #
##################################
##### Create triplets (2.5th percentile, point estimate, 97.5th percentile)
### Step 1a: fit model
logmod <- glm(y ~ x1 + x2, 
              data=df, 
              family="binomial")

### Step 1b: extract estimates for tables, betas, variance-covariance matrix
summary(logmod)                 # Betas for table
exp(logmod$coefficients)        # ORs for table
exp(confint(logmod, level=.95)) # OR CIs for table
pes <- coef(logmod)             
vcov <- vcov(logmod)  

### Step 2a: set up counterfactuals (x values at which to probe model)
x1range <- seq(-3, 3, by=.1) # range of hypothetical x1 values
scenlen <- length(x1range)   # number of x1 scenarios per level of x2
xi1 <- xi2 <- as.data.frame(matrix(NA,
                            nrow=scenlen,
                            ncol=length(pes))) #xi1 is data for x2=low, xi2 is data for x2=high
names(xi1) <- names(xi2) <- names(pes)
xi1$`(Intercept)` <- xi2$`(Intercept)` <- rep(1,scenlen)
xi1$x1 <- xi2$x1 <- x1range
xi1$x2 <- rep(-1, scenlen)
xi2$x2 <- rep(1, scenlen)
matxi1 <- as.matrix(xi1)
matxi2 <- as.matrix(xi2)

### Step 2b: generate predictions from the model we fit from our data
pred1 <- predict(logmod, newdata=xi1, type="response")
pred2 <- predict(logmod, newdata=xi2, type="response")

### Step 3a: simulate a range of plausible models based on distribution of parameters from model fit 
sims <- 10000
simparams <- mvrnorm(sims,pes,vcov)

### Step 3b: create triplets for each set of scenarios
lo <- pe <- hi <- rep(NA,scenlen)
# for x2=low
trip1 <- data.frame(lo=lo, pe=pred1, hi=hi) # create triplets for x2=low
for(i in 1:nrow(matxi1)){ 
  simmu <- simparams%*%matxi1[i,]  # result of XBeta for all simulated models
  simy <- 1/(1 + exp(-simmu))      # logistic function - convert result to a probability
#  simy <- exp(simmu)              # exponential function - inverse of log-link, uncomment if Poisson
  simy <- sort(simy)
  
  trip1$lo[i] <- quantile(simy,.025)
  trip1$hi[i] <- quantile(simy,.975)
}
trip2 <- data.frame(lo=lo, pe=pred2, hi=hi) # create triplets for x2=high
for(i in 1:nrow(matxi2)){ 
  simmu <- simparams%*%matxi2[i,] 
  simy <- 1/(1 + exp(-simmu)) 
  simy <- sort(simy)
  
  trip2$lo[i] <- quantile(simy,.025)
  trip2$hi[i] <- quantile(simy,.975)
}
dfplot1 <- data.frame(x1=xi1$x1, x2=xi1$x2, 
                      pe=trip1$pe, lo=trip1$lo, hi=trip1$hi) # attach triplets to predictors
dfplot2 <- data.frame(x1=xi2$x1, x2=xi2$x2, 
                      pe=trip2$pe, lo=trip2$lo, hi=trip2$hi)
dfplot <- rbind(dfplot1, dfplot2) # create a single plotting dataset
dfplot$x2 <- as.factor(dfplot$x2)
levels(dfplot$x2) <- c("High (+1 SD)","Low (-1 SD)") 

##############################
# Create graphics and tables #
##############################
### Plot triplets as fit line plus region of uncertainty
(ex1p <- ggplot(data=dfplot, mapping=aes(x=x1, y=pe, linetype=x2)) +
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=.10) +
    geom_line() +
    labs(linetype=expression(x[2])) +
    xlab(expression(x[1])) +
    ylab("Expected Probability") +
    scale_y_continuous(breaks=seq(0, 1, .2), minor_breaks=seq(0, 1, .1), limits=c(0, 1)) +
    theme_bw()
)
h <- 6
ggsave(file="Images/p_by_x1x2.png", 
       height=h,
       width=h*1.619)

### Calculate predicted values for tables
dftab <- data.frame(x0=rep(1, times=9),
                    x1=rep(c(-1, 0, 1), each=3), 
                    x2=rep(c(-1, 0, 1), times=3))
dftab$yev <- predict(logmod, newdata=dftab, type="response") # point estimates
dftab

### Calculate 95% confidence intervals for individual point estimates
names(dftab)[1] <-"(Intercept)"
simmu <- simparams%*%t(as.matrix(dftab)[,1:3]) 
simy <- 1/(1 + exp(-simmu))

dftab$lo <- apply(simy, 2, quantile, .025) 
dftab$hi <- apply(simy, 2, quantile, .975)

### Have R output these values
paste0("x1 = ", dftab$x1,
       "; x2 = ", dftab$x2,
       "; Predicted probability: ", round(dftab$yev, 2), 
       " [",round(dftab$lo, 2), ", ", round(dftab$hi, 2), "]")