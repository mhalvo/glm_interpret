#########################################
# sample_code_emp_ex.R                  #
#    Code used to generate figures      #
#     and table estimates for empirical #
#     paper example                     #
#########################################
# Max Halvorson #
#   2020-05-07  #
#################

#######################
# Setup and data load #
#######################
rm(list=ls())
library('MASS')
library('ggplot2')
library('psych')
options(max.print=2000, scipen=7, tibble.print_min=70)

setwd('C:/Users/Max Halvorson/Dropbox/GLM Main Effects Paper')

dat <- read.table("https://raw.githubusercontent.com/mhalvo/glm_interpret/master/ex_clean.txt", header = T)
dat$gen <- dat$gen - 1
dat$purg_new <- dat$purg - 1


##################################
# Analysis for graphs and tables #
##################################
### Create triplets (2.5th percentile, point estimate, 97.5th percentile)
# Step 1a: fit model
pmod <- glm(negpysum ~ as.factor(gen) + purg_new + pre + ss, 
            family=poisson(link="log"),
            data = dat)
summary(pmod)
# Step 1b: extract estimates for tables, betas, variance-covariance matrix
exp(coef(pmod))                # ORs, IRRs for table
exp(confint(pmod, level=.95))  # OR, IRR CIs for table
p_pes <- coef(pmod)
p_vcov <- vcov(pmod)

# Step 2a: set up counterfactuals (x values at which to probe model)
purgrange <- seq(0, 3, by=.25) # hypothetical purg levels
scenlen <- length(purgrange)   # number of scenarios per levels of other predictors
desc <- psych::describe(dat)
gen_m <- 0
gen_f <- 1
pre_m <- desc$mean[7]
ss_m <- desc$mean[9]
# male
df_male <- data.frame(int = 1,
                      gen = rep(gen_m, scenlen),
                      purg_new = purgrange,
                      pre = rep(pre_m, scenlen),
                      ss = rep(ss_m, scenlen))
# female
df_female <- data.frame(int = 1,
                        gen = rep(gen_f, scenlen),
                        purg_new = purgrange,
                        pre = rep(pre_m, scenlen),
                        ss = rep(ss_m, scenlen))
names(df_male)[1] <- names(df_female)[1] <- "(Intercept)"

# Step 2b: generate expected values from our data model 
pred_male <- predict(pmod, newdata=df_male, type="response")
pred_female <- predict(pmod, newdata=df_female, type="response")

### Step 3a: simulate a range of plausible models based on distribution of parameters from model fit 
sims <- 10000
simparams <- mvrnorm(sims,p_pes,p_vcov)

# Step 3b: create triplets for each set of scenarios
matdfm <- as.matrix(df_male)
matdff <- as.matrix(df_female)
lo <- pe <- hi <- rep(NA,scenlen)

tripm <- data.frame(lo=lo, pe=pred_male, hi=hi) # create triplets for male
for(i in 1:scenlen){
  simmu <- simparams %*% matdfm[i,]    # result of Xbeta for count
  simy <- exp(simmu)                  # inverse of log-link - convert to a count
  
  tripm$lo[i] <- quantile(simy,.025)
  tripm$hi[i] <- quantile(simy,.975)
}

tripf <- data.frame(lo=lo, pe=pred_female, hi=hi) # create triplets for female
for(i in 1:scenlen){ 
  simmu <- simparams %*% matdff[i,]  
  simy <- exp(simmu)
  
  tripf$lo[i] <- quantile(simy,.025)
  tripf$hi[i] <- quantile(simy,.975)
}

dfplotm <- cbind(matdfm, tripm)
dfplotf <- cbind(matdff, tripf)
dfplot <- rbind(dfplotm, dfplotf)
dfplot$gen <- as.factor(dfplot$gen)
levels(dfplot$gen) <- c("Male","Female")

##############################
# Create graphics and tables #
##############################
### Plot triplets as fit line plus region of uncertainty
(ex2p <- ggplot(data=dfplot, mapping=aes(x=purg_new, y=pe, linetype=gen)) +
    geom_line() +
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=.10) +
    labs(linetype="Gender") +
    xlab("Positive Urgency") +
    ylab("Expected Count of Alcohol Consequences") +
    ylim(0, 17) +
    scale_x_continuous(breaks=c(0,1,2,3),
                       labels=c("1 - Disagree\n Strongly","2 - Disagree\n Somewhat",
                                "3 - Agree\n Somewhat","4 - Agree\n Strongly")) +
    theme_bw()
)  
h <- 6
ggsave("Images/alcons_by_gender.png",
       height=h, 
       width=h*1.619)

### Calculate predicted values for tables
dftab <- data.frame(int=rep(1,8),
                    gen=rep(c(0, 1), times=4),
                    purg_new=rep(0:3, each=2),
                    pre=pre_m,
                    ss=ss_m)
dftab$yev <- predict(pmod, newdata=dftab, type="response") # point estimates
dftab

### Calculate 95% confidence intervals for individual point estimates
names(dftab)[1] <-"(Intercept)"
dftabm <- as.matrix(dftab[,-6])
simmu <- simparams%*%t(dftabm) 
simy <- exp(simmu)

dftab$lo <- apply(simy, 2, quantile, .025) 
dftab$hi <- apply(simy, 2, quantile, .975)

### Have R output these values
paste0("Gender = ", dftab$gen,
       "; Pos Urg = ", round(dftab$purg_new, 2),
       "; Planning = ", round(dftab$pre, 2),
       "; Sen Seek = ", round(dftab$ss, 2),
       "; Predicted count: ", round(dftab$yev, 2), 
       " [",round(dftab$lo, 2), ", ", round(dftab$hi, 2), "]")