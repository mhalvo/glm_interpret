##########################################
# sample_code_generic.R                  #
#    Code used to generate figures       #
#     and table estimates for simulated  #
#     paper example: meant to be adapted #
#     for individual use cases           #
##########################################
# Max Halvorson #
#   2020-03-06  #
#################

rm(list=ls()) 
library('MASS')       # for mvrnorm
library('tidyverse')  # for ggplot
library('psych')      # for describe
options(max.print=2000, scipen=7, tibble.print_min=70)

setwd('C:/Users/Max Halvorson/Documents/GitHub/glm_main_effects')

##########################
# Simulated data example #
##########################
### Simulate data
set.seed(1337) 
b0 <- 1        # Intercept
b1 <- -1       # X1 Effect
b2 <- 1        # X2 Effect
n <- 500       # Sample Size
mu <- rep(0,2) # Specify means
S <- matrix(c(  1,-.1,
              -.1,  1),nrow=2,ncol=2) # Specify covariance matrix
sigma <- 1                            # Level 1 error
rawvars <- mvrnorm(n=n, mu=mu, Sigma=S) # Simulate predictors from a multivariate normal distribution

id <- seq(1:n)
eij <- rep(rnorm(id, 0, sigma)) # Add normally distributed error
xb <- (b0) + (b1)*rawvars[,1] + (b2)*rawvars[,2] + eij

pr <- 1 / (1 + exp(-xb))        # pass through logistic function to convert xb to probabilities
y <- rbinom(n, 1, pr)       # bernoulli response variable
df <- data.frame(id=id,
                 y=y,
                 x1=rawvars[,1],
                 x2=rawvars[,2])

##################################
# Analysis for graphs and tables #
##################################
##### Create triplets (2.5th percentile, point estimate, 97.5th percentile)
### Step 1: fit model
logmod <- glm(y ~ x1 + x2, 
              data=df, 
              family="binomial")

### Step 2: extract estimates for tables, betas, variance-covariance matrix
summary(logmod)                 # Betas for table
exp(logmod$coefficients)        # ORs for table
exp(confint(logmod, level=.95)) # OR CIs for table
pes <- coef(logmod)             
vcov <- vcov(logmod)  

### Step 3: set up counterfactuals (x values at which to probe model)
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

### Step 4: simulate many plausible models based on distribution of parameters from model fit in Step 1 
sims <- 10000
simparams <- mvrnorm(sims,pes,vcov)

### Step 5: create triplets for each set of scenarios
lo <- pe <- hi <- rep(NA,scenlen)
# for x2=low
trip1 <- data.frame(lo=lo, pe=pe, hi=hi) # create triplets for x2=low
for(i in 1:nrow(matxi1)){ 
  simmu <- simparams%*%matxi1[i,]  # result of XBeta for all simulated models
  simy <- 1/(1 + exp(-simmu))      # logistic function - convert result to a probability
  simy <- sort(simy)
  
  trip1$lo[i] <- quantile(simy,.025)
  trip1$pe[i] <- quantile(simy,.50)
  trip1$hi[i] <- quantile(simy,.975)
}
trip2 <- data.frame(lo=lo, pe=pe, hi=hi) # create triplets for x2=high
for(i in 1:nrow(matxi2)){ 
  simmu <- simparams%*%matxi2[i,] 
  simy <- 1/(1 + exp(-simmu)) 
  simy <- sort(simy)
  
  trip2$lo[i] <- quantile(simy,.025)
  trip2$pe[i] <- quantile(simy,.50)
  trip2$hi[i] <- quantile(simy,.975)
}
dfplot1 <- data.frame(x1=xi1$x1, x2=xi1$x2, 
                      pe=trip1$pe, lo=trip1$lo, hi=trip1$hi) # attach triplets to predictors

dfplot2 <- data.frame(x1=xi2$x1, x2=xi2$x2, 
                      pe=trip2$pe, lo=trip2$lo, hi=trip2$hi)
dfplot <- rbind(dfplot1, dfplot2) # create a single plotting dataset
dfplot$x2 <- as.factor(dfplot$x2)
levels(dfplot$x2) <- c("High (+1 SD)","Low (-1 SD)") 

##### Create graphics and  tables
### Step 6: plot triplets as fit line plus region of uncertainty
(ex1p <- ggplot(data=dfplot, mapping=aes(x=x1, y=pe, linetype=x2)) +
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=.10) +
    geom_line() +
    labs(linetype="X2") +
    xlab("X1") +
    ylab("Expected Probability") +
    scale_y_continuous(breaks=seq(0, 1, .2), minor_breaks=seq(0, 1, .1), limits=c(0, 1)) +
    theme_bw()
)
h <- 6
ggsave(file="p_by_x1x2.png",
       height=h,
       width=h*1.619)

### Step 7: calculate values for tables
dftab <- data.frame(x1=rep(c(-1, 0, 1), each=3), 
                    x2=rep(c(-1, 0, 1), times=3))
dftab$yev <- predict(logmod, newdata=dftab, type="response") # plug-and-chug new predictors through model
dftab
