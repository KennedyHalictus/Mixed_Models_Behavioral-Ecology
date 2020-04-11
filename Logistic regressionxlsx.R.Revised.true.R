## try to keep all the libraries on the top
## it's easier to keep track of what you're doing

library(readxl)
library(zoo) 
library(dplyr)

library(tidyr)
library(lme4) #?GLMM
library(gridExtra)
library(ggplot2); theme_set(theme_bw())
getwd()
paireddata <- read_xlsx("spring bees.paired aurata.xlsx")

paireddata$nest <- na.locf(paireddata$nest)
paireddata$day <- na.locf(paireddata$day)

paireddata2 <- paireddata %>% 
  mutate(
    nest=na.locf(nest),
    day=na.locf(day) 
  ) %>%
  select(-time) %>%
  gather(key, value, -nest, -day, -beeid, -observation) %>% 
  mutate(
    value=ifelse(is.na(value), 0, value)
  )

#Explaratory analyses

beesummed <- paireddata2 %>%
  group_by(day, key) %>%  #grouping the data by day 1-12. Summing variables instead of days?
  summarize(
    prop=mean(value)  #proportion? #yes #per nest per day???
  )

#Logistic Regression

library(ggplot2)

g0 <- ggplot(beesummed) +
  geom_line(aes(day, prop)) + 
  geom_point(aes(day, prop)) + 
  ylab("Mean value of observed behaviours") + 
  facet_wrap(~key, nrow=6) +
  theme(
    strip.background = element_blank()
  ) + xlim(0, 12.5)
g0

paireddata3 <- paireddata2 %>%
  filter(key !="gc")  #no activity on gc.Its removal from analysis. Preliminary aanalysis without random effects on the nests

#SEPARATING FITS FOR EACH VARIABLE

behavior <-unique(paireddata3$key)
glmfitlist <-vector('list', length(behavior)) 
names(glmfitlist) <- behavior
  
for (i in 1:length(behavior)) {
  tmpdata <-paireddata3 %>%
    filter(key==behavior[i]) %>%
    mutate(day=day-1)
    
  gfit <- try(glm(value~1 + day, family="binomial", data=tmpdata))
    
  if (!inherits(gfit, "try-error")) {
    glmfitlist[[i]] <- gfit
  }
}
  
#mixed models.Estimates the effects of one or more explanatory variables on a response variable
#Comparing parameter estimates between behaviors
glmest <-                            #glmestimate
  lapply(glmfitlist, function(x) {
    
    cc <- confint(profile(x))
    
    data.frame(
      par=c("intercept", "slope"),
      est=coef(x),
      lwr=cc[,1],
      upr=cc[,2]
    )
    }) %>%
    bind_rows(.id="behavior")

glmest1 <- glmest %>%
  filter(par=="intercept") %>%
  arrange(est) %>%
  mutate(behavior=factor(behavior, level=behavior)) #Intro behavior as a factor.
glmest2 <- glmest %>%
  filter(par=="slope") %>%
  arrange(est) %>%
  mutate(behavior=factor(behavior, level=behavior))

g01 <- ggplot(glmest1) +                      #intercept
  geom_point(aes(behavior, est)) +
  geom_errorbar(aes(behavior, ymin=lwr, ymax=upr)) +          
  geom_hline(yintercept=0, lty=2) +
  ylab("baseline log-odds at day 1") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

g01

g02 <- (g01 %+% glmest2) +                 #Intercept + Slope
  ylab("log-odds ratio (per day)")


grid.arrange(g01, g02, nrow=1)

#Imposing mixed models 

glmerfitlist <- vector('list', length(behavior)) 
names(glmerfitlist) <- behavior
  
for (i in 1:length(behavior)) {        
  tmpdata <- paireddata3 %>%
    filter(key==behavior[i]) %>%
    mutate(day=day-1)
    
  gfit <- try(glmer(value~1 + day +(day|nest), family="binomial", data=tmpdata))
    
  if (!inherits(gfit, "try-error")) {
    glmerfitlist[[i]] <- gfit
  }
}    
#Wald confidence intervals(for final analysis, use profile or bootstrap)
  
glmerest <-
  lapply(glmerfitlist,function(x) {
    #again, doing this silences the the confidenceinterval
    cc <- confint(x, parm=c("(Intercept)", "day"), method="Wald")
    
    data.frame(
      par=c("intercept", "slope"), #
      est=fixef(x),
      lwr=cc[,1],
      upr=cc[,2]
      
    )
  }) %>%
  bind_rows(.id="behavior")
  
glmerest1 <- glmerest %>%
  filter(par=="intercept") %>%
  arrange(est) %>%
  mutate(behavior=factor(behavior, level=behavior))

glmerest2 <- glmerest %>%        
  filter(par=="slope") %>%
  arrange(est) %>%
  mutate(behavior=factor(behavior, level=behavior))

glmcomb1 <- mutate(glmerest1, fit="glmer") %>%     #intercepts combined
  bind_rows(mutate(glmest1, fit="glm"))

glmcomb2 <- mutate(glmerest2, fit="glmer") %>%    #slopes combined
  bind_rows(mutate(glmest2, fit="glm"))


g03 <- ggplot(glmcomb1) +        #ggplot intercepts combined
  geom_point(aes(behavior, est, col=fit), position = position_dodge(0.5)) +
  geom_errorbar(aes(behavior, ymin=lwr, ymax=upr, col=fit), width=0, position = position_dodge(0.5)) +
  geom_hline(yintercept=0, lty=2) +
  ylab("baseline log-odds at day 1") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.title = element_blank(),
    legend.position = "top"
    )
g03

g04 <- (g03 %+% glmcomb2) +
  ylab("log-odds ratio(per day)")

grid.arrange(g03, g04, nrow=1)

#ESTIMATED PROBABILITY VS WHAT WE OBSERVED
predlist <- vector('list',length(behavior))
names(predlist) <- behavior

for (i in 1:length(behavior)) {
  gg <- glmerfitlist[[i]]
  
  predlist[[i]] <- data.frame(
    day=1:13,
    pred=predict(gg, newdata=data.frame(day=0:12), re.form=NA, type="response")
  )
}

preddata <- predlist %>%
  bind_rows(.id="key")

## should include confidence intervals on the estimated probabilities for the final analysis...
g0 +
  geom_line(data=preddata, aes(day, pred), col="red")

## you could also consider a model like this
## and model multiple behaviors at the same time
## and try to measure the underlying correlation among different behaviours
## this is one of the simplest models of this kind you could consider
## almost better to go Bayesian if you want to further pursue this direction
## a lot of potentially interesting ways you can go with this (and probably more realistic)...
## showing you an example with two behaviours (fe and ab)

gfit <- glmer(value~-1 + key + day:key +
                (-1 + key|nest) +
                (-1 + day:key|nest), data=filter(paireddata3, key%in%c("an", "ab")),
                family=binomial)

VarCorr(gfit)

## ignore NAs
confint(gfit, method="Wald", oldNames=FALSE)

summary(gfit)
