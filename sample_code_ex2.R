#########################################
# sample_code_ex2.R                     #
#    Code used to generate figures      #
#     and table estimates for empirical #
#     paper example                     #
#########################################
# Max Halvorson #
#   2020-03-06  #
#################

#######################
# Setup and data load #
#######################
rm(list=ls())
library('MASS')       # for mvrnorm
library('tidyverse')  # for ggplot
library('psych')      # for describe
library('pscl')       # for zeroinfl
options(max.print=2000, scipen=7, tibble.print_min=70)

dat <- read.table("https://raw.githubusercontent.com/mhalvo/glm_main_effects/master/ex2_clean.txt", header = T)
dat$gen <- dat$gen - 1
dat$purg_new <- dat$purg - 1

##################################
# Analysis for graphs and tables #
##################################

### Create triplets
# Step 1: fit model
zimod <- zeroinfl(negpysum ~ as.factor(gen) + purg_new + pre | ss, data = dat)
summary(zimod)

# Step 2: extract betas, vcov
exp(coef(zimod))                # ORs, IRRs for table
exp(confint(zimod, level=.95))  # OR, IRR CIs for table
zi_pes <- coef(zimod)
zi_vcov <- vcov(zimod)

# Step 3: set up counterfactuals
purgrange <- seq(0, 3, by=.25) # hypothetical purg levels
scenlen <- length(purgrange)   # number of scenarios per levels of other predictors
desc <- psych::describe(dat)
gen_m <- 0
gen_f <- 1
pre_m <- desc["pre","mean"]
ss_m <- desc["ss", "mean"]
# male
df_male_c <- data.frame(int_c = 1,
                        gen = rep(gen_m, scenlen),
                        purg_new = purgrange,
                        pre = rep(pre_m, scenlen))
df_male_z <- data.frame(int_z = 1,
                        ss = rep(ss_m, scenlen))
# female
df_female_c <- data.frame(int_c = 1,
                          gen = rep(gen_f, scenlen),
                          purg_new = purgrange,
                          pre = rep(pre_m, scenlen))
df_female_z <- data.frame(int_z = 1,
                          ss = rep(ss_m, scenlen))

# Step 4: simulate plausible models from model fit in Step 1 
sims <- 10000
simparams <- mvrnorm(sims,zi_pes,zi_vcov)
simparamsc <- simparams[,1:4]
simparamsz <- simparams[,5:6]

matdfmc <- as.matrix(df_male_c)
matdfmz <- as.matrix(df_male_z)

matdffc <- as.matrix(df_female_c)
matdffz <- as.matrix(df_female_z)

# Step 5: create triplets for each set of scenarios
lo <- pe <- hi <- rep(NA,scenlen)
tripm <- data.frame(lo=lo, pe=pe, hi=hi) # create triplets for male
for(i in 1:scenlen){ 
  # count portion
  simmuc <- simparamsc %*% matdfmc[i,]  # result of Xbeta for count
  simyc <- exp(simmuc)                  # inverse of log-link - convert to a count
  # zero portion
  simmuz <- simparamsz %*% matdfmz[i,]  # result of Xbeta for zero
  simyz <- 1/(1 + exp(-simmuz))         # inverse of logit-link - convert to a P(0)
  # integrate
  simy <- simyc * (1 - simyz) + 0 * simyz # mixture model: weighted average of count and zero expected values
  
  tripm$lo[i] <- quantile(simy,.025)
  tripm$pe[i] <- quantile(simy,.50)
  tripm$hi[i] <- quantile(simy,.975)
}
tripf <- data.frame(lo=lo, pe=pe, hi=hi) # create triplets for female
for(i in 1:scenlen){ 
  # count portion
  simmuc <- simparamsc %*% matdffc[i,]  
  simyc <- exp(simmuc)                  
  # zero portion
  simmuz <- simparamsz %*% matdffz[i,]  
  simyz <- 1/(1 + exp(-simmuz))         
  # integrate
  simy <- (1 - simyz)*simyc + (simyz)*0
  
  tripf$lo[i] <- quantile(simy,.025)
  tripf$pe[i] <- quantile(simy,.50)
  tripf$hi[i] <- quantile(simy,.975)
}

dfplotm <- cbind(matdfmc, matdfmz, tripm)
dfplotf <- cbind(matdffc, matdffz, tripf)
dfplot <- rbind(dfplotm, dfplotf)
dfplot$gen <- as.factor(dfplot$gen)
levels(dfplot$gen) <- c("Male","Female")

### Step 6: plot triplets
(ex2p <- ggplot(data=dfplot, mapping=aes(x=purg_new, y=pe, linetype=gen)) +
    geom_line() +
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=.10) +
    labs(linetype="Gender") +
    xlab("Positive Urgency") +
    ylab("Expected Count of Alcohol Consequences") +
    ylim(0, 16) +
    scale_x_continuous(breaks=c(0,1,2,3),
                       labels=c("1 - Disagree\n Strongly","2 - Disagree\n Somewhat",
                                "3 - Agree\n Somewhat","4 - Agree\n Strongly")) +
    theme_bw()
)  

h <- 6
ggsave("alcons_by_purggen.png",
       height=h, 
       width=h*1.619)

### Step 7: calculate values for tables
dftab <- data.frame(gen=rep(c(0, 1), times=4),
                    purg_new=rep(0:3, each=2),
                    pre=pre_m,
                    ss=ss_m)
dftab$yev <- predict(zimod, newdata=dftab, type="response")
dftab
