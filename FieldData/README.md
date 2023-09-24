## Behaviour Data
```{r Load data}
library(tidyverse)
####Set working directory 
setwd("~/Documents/Trancriptomics A. octoarticulatus/FieldData")
##data <- read_csv("~/Documents/Trancriptomics A. octoarticulatus/FieldData/FieldDataAllomerus.csv")
View(data)

#Fix factors, etc.
data$plant <- as.factor(data$plant)
data$gh.species <- as.factor(data$gh.species)
data$gh.age <- as.factor(data$gh.age)

library(hms)
#Fix time
data$time <- hms(data$time)
data$time2 <- as.numeric(data$time)

###Explore data
hist(sqrt(ao_data$activity.avg))
```

```{r colony consistency}
install.packages("maditr")
install.packages("MASS")
install.packages("pscl")
library(maditr)
library(MASS)
library(pscl)

#Reorganize to model data at two time points
data_wide <- dcast(data, plant~trial, value.var="activity.avg", drop=FALSE)

p1 <- ggplot(data=data_wide, aes(x=sqrt(`1`), y=sqrt(`3`)))+
      geom_point()+
      geom_smooth(method="lm")
p1

p2 <- ggplot(data=data_wide, aes(x=sqrt(`1`), y=sqrt(`2`)))+
      geom_point()+
      geom_smooth(method="lm")
p2

lm1 <- lm(sqrt(`1`) ~ sqrt(`3`), data=data_wide)
summary(lm1) #Same as what Chris gets

nb_model <- glm.nb(round(`1`)~round(`3`), data=data_wide)
summary(nb_model) #NS

zero_inflated_model_poisson <- zeroinfl(round(`1`) ~ round(`3`)|1, dist="poisson", data=data_wide)
summary(zero_inflated_model_poisson) #NS

zero_inflated_model_nb <- zeroinfl(round(`1`) ~ round(`3`)|1, dist="negbin", data=data_wide)
summary(zero_inflated_model_nb) #NS

poisson_model <- glm(round(`1`)~round(`3`), data=data_wide, family="poisson")
summary(poisson_model)
plot(poisson_model)

vuong(nb_model, zero_inflated_model_nb) #Zero-inflated model not better
vuong(poisson_model, zero_inflated_model_poisson) #Zero-inflated Poisson model is better

```

```{r make data long}
data_long <- gather(data, minute, ants, min1:min5, factor_key=TRUE)
data_long$trial_minute <- paste0(data_long$trial, "_", data_long$minute)
data_long$plant <- as.factor(data_long$plant)
hist(sqrt(data_long$ants))

#Model trial 1 data with all variables
install.packages("lme4")
library(lme4)

#This models the square-root transformed average number of ant bodyguards
lm1 <- lm(sqrt(activity.avg)~temp.c +time2 + dom.no + gh.length.mm + gh.age + gh.species, data=subset(data, trial == "1"))
summary(lm1)
plot(lm1)

lm2 <-  lm(sqrt(activity.avg)~time2 + gh.species, data=subset(data, trial == "1"))
summary(lm2)
plot(lm2)

#This models the count at each minute, with plant ID as a random effect (gets the same answers)
lmm1 <- lmer(sqrt(ants)~trial_minute + temp.c +time2 + dom.no + gh.length.mm + gh.age + gh.species + (1|plant), data=subset(data_long, trial == "1"))
summary(lmm1)
plot(lmm1)

#Drop non-significant variables, and re-fit
lmm2 <- lmer(sqrt(ants)~trial_minute + time2  + gh.species + (1|plant), data=subset(ao_data_long, trial == "1"))
summary(lmm2)
plot(lmm2)

#Model trial 3 data with all variables, then drop non-significant variables and re-fit
lm3 <-  lm(sqrt(activity.avg)~temp.c + dom.no + dom.size.mm+gh.length.mm+time2 + gh.species, data=subset(data, trial == "3"))
summary(lm3)
plot(lm3)

lm4 <-  lm(sqrt(activity.avg)~time2, data=subset(data, trial == "3"))
summary(lm4)
plot(lm4)

#Visualize time-of-day effects
p3 <- ggplot(data=data, aes(x=time, y=sqrt(activity.avg)))+
      geom_point()+
      scale_x_time()+
      facet_wrap(~trial)
p3

#Visualize grasshopper species effect (in trial 1 only)
p4 <- ggplot(data=subset(data, trial == "1"), aes(x=gh.species, y=sqrt(activity.avg)))+
      geom_boxplot()
p4
```
