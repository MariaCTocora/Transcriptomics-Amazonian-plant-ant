### Field data analysis

This markdown document analyzes the field bioassay and herbivory data from:

MC Tocora, C Reid, H Xue, ME Frederickson. Going viral: transcriptomics of an Amazonian ant-plant symbiosis reveals viral infections correlate with ant ‘bodyguard’ behavior. In preparation.

First, let's load packages. 

```{r load packages, warning=FALSE, message=FALSE}
library(tidyverse)
library(lubridate)
library(cowplot)
library(lme4)
library(png)
library(car)
library(magick)
```

Next, let's read in the field data.

```{r load and clean data, warning=FALSE, message=FALSE}
#Read in data
df <- read.csv("FieldData.csv")

#Label colonies that were categorized as high-quality or low-quality bodyguards in RNA-Seq project
df$bodyguard.cat <- ifelse(df$plant %in% c(103, 199, 28, 207, 187, 18, 109), "Low", ifelse(df$plant %in% c(107, 145, 136, 201, 193, 127, 125), "High", NA))

#Make bodyguarding a factor
df$bodyguard.cat <- as.factor(df$bodyguard.cat)
```

Choose the colors used in the figures. 

```{r set figure colors, warning=FALSE, message=FALSE}
fig_colors <- c("blue", "deeppink", "black")
```

Analyze ant bodyguarding behavior. 

```{r ant bodyguard behavior, warning=FALSE, message=FALSE}
#Make wide data long
df_long <- gather(df, bioassay_minute, ants, min1:min5)
#Make the minutes elapsed in the bioassay numeric
df_long$bioassay_minute_numeric <- as.numeric(gsub(".*?([0-9]+).*", "\\1", df_long$bioassay_minute))

p1 <- ggplot()+geom_line(data=df_long, aes(x=bioassay_minute_numeric, y=ants, group=plant), alpha=0.5, color="grey")+geom_line(data=subset(df_long, !is.na(bodyguard.cat)), aes(x=bioassay_minute_numeric, y=ants, color=bodyguard.cat, group=plant))+theme_cowplot()+xlab("Time (min.)")+ylab("Ant bodyguards (no.)")+labs(color="Bodyguard activity")+scale_color_manual(values=fig_colors)+theme(legend.position=c(0.1, 0.85))+scale_y_continuous(limits=c(0,50))
p1

#Does bodyguarding activity differ between the two groups? 
lm1 <- lm(activity.avg~bodyguard.cat, data=subset(df, !is.na(bodyguard.cat)))
summary(lm1)
Anova(lm1, type=3)
min.ants.high <- min(subset(df, bodyguard.cat == "High")$activity.avg)
max.ants.low <- max(subset(df, bodyguard.cat == "Low")$activity.avg)

#Plot(lm1)
p1 + annotate(geom = "text", x=4.5, y=48, label=paste0("p = ", round(summary(lm1)$coefficients[2,4], 3)))

```

Analyze herbivory.

```{r herbivory, warning=FALSE, message=FALSE}
#Does ant behavior in the bioassay predict herbivory? 
lm2 <- lm(herbivory.rank~activity.avg, data=df)
summary(lm2)
Anova(lm2, type=3)
#plot(lm2)

p2 <- ggplot(data=df)+geom_point(aes(x=activity.avg, y=herbivory.rank, color=bodyguard.cat), size=2)+theme_cowplot()+geom_smooth(aes(x=activity.avg, y=herbivory.rank),method="lm", se=FALSE, color="black")+scale_color_manual(values=fig_colors)+xlab("Ant bodyguards (mean no.)")+ylab("Herbivory (rank)")+theme(legend.position="none")+annotate(geom = "text", x=28, y=55, label=paste0("p = ", round(summary(lm2)$coefficients[2,4], 3)))
p2

#Does herbivory differ between the two groups? 
lm3 <- lm(herbivory.rank~bodyguard.cat, data=df)
summary(lm3)
Anova(lm3, type=3)

p3 <- ggplot(data=subset(df, !is.na(bodyguard.cat)))+geom_boxplot(aes(x=bodyguard.cat, y=herbivory.rank, fill=bodyguard.cat))+theme_cowplot()+scale_fill_manual(values=fig_colors)+xlab("Bodyguard activity")+ylab("Herbivory (rank)")+theme(legend.position="none")+annotate(geom = "text", x=1.5, y=55, label=paste0("p = ", round(summary(lm3)$coefficients[2,4], 3)))
p3
```

Make multi-panel plot. 

```{r multipanel plot, warning=FALSE, message=FALSE}
#Read in png figure
image <- readPNG("Figure_1.png")

#Set png figure up in ggplot
system_fig <- ggdraw() +
  draw_image(
    image, scale = 0.85, x = 1, y = 0.1,
    hjust = 1, halign = 1, valign = 0
  )  

#Assemble multi-panel figure
p4 <- plot_grid(system_fig, p1, p2, p3, labels="AUTO")
p4

#Save plot
save_plot("Figure1.pdf", p4, base_width=8, base_height=6)
```

Analyze other predictors of ant bodyguarding behavior.

```{r other predictors, warning=FALSE, message=FALSE}
#Make time a time 
df$time.fixed <- lubridate::hms(df$time)
df$time.length.min <- time_length(df$time.fixed - hms("00:00:00"), unit = "minute")

#Model the data
lm4 <- lm(activity.avg~time.length.min+temp.c+dom.no+gh.mass.g+gh.species+gh.age, data=df)
summary(lm4)
Anova(lm4, type=3)

p5 <- ggplot(data=df)+geom_point(aes(x=time.fixed, y=activity.avg, color=bodyguard.cat))+geom_smooth(aes(x=time.fixed, y=activity.avg), method="lm", se=FALSE)+theme_cowplot()+scale_color_manual(values=fig_colors)+xlab("Time of day")+ylab("Ant bodyguards (no.)")+theme(legend.position="none")+scale_x_time(labels=scales::time_format("%H:%M"), limits=c(hms("8:00:00"), hms("18:00:00")))+annotate(geom = "text", x=hms("13:00:00"), y=30, label=paste0("p = ", round(Anova(lm4)[[4]][[1]], 3)))
p5

p6 <- ggplot(data=df)+geom_point(aes(x=temp.c, y=activity.avg, color=bodyguard.cat))+theme_cowplot()+scale_color_manual(values=fig_colors)+xlab(expression("Temperature " (degree*C)))+ylab("Ant bodyguards (no.)")+theme(legend.position="none")+annotate("text", x=25, y=30, label=paste0("p = ", round(Anova(lm4)[[4]][[2]], 3)))+scale_x_continuous(limits=c(23,27))
p6

p7 <- ggplot(data=df)+geom_point(aes(x=dom.no, y=activity.avg, color=bodyguard.cat))+theme_cowplot()+scale_color_manual(values=fig_colors)+xlab("Domatia (no.)")+ylab("Ant bodyguards (no.)")+theme(legend.position="none")+annotate(geom = "text", x=62.5, y=30, label=paste0("p = ", round(Anova(lm4)[[4]][[3]], 3)))+scale_x_continuous(limits=c(0,125))
p7

p8 <- ggplot(data=df)+geom_point(aes(x=gh.mass.g, y=activity.avg, color=bodyguard.cat))+theme_cowplot()+scale_color_manual(values=fig_colors)+xlab("Grasshopper size (g)")+ylab("Ant bodyguards (no.)")+theme(legend.position="none")+annotate(geom = "text", x=0.025, y=30, label=paste0("p = ", round(Anova(lm4)[[4]][[4]], 3)))+scale_x_continuous(limits=c(0,0.05))
p8

p9 <- ggplot(data=subset(df, gh.age == "adult" | gh.age == "juvenile"))+geom_boxplot(aes(x=gh.age, y=activity.avg))+theme_cowplot()+scale_color_manual(values=fig_colors)+xlab("Grasshopper age")+ylab("Ant bodyguards (no.)")+theme(legend.position="none")+scale_x_discrete(labels=c("Adult", "Juvenile"))+annotate(geom = "text", x=1.5, y=30, label=paste0("p = ", round(Anova(lm4)[[4]][[6]], 3)))
p9

#Trim white space from grasshopper names
df$gh.species <- trimws(df$gh.species)

p10 <- ggplot(data=subset(df, gh.species == "Orphulella concinnula" | gh.species == "Orphulella fluvialis" | gh.species == "Orphulella punctata" | gh.species == "Orphulella sp."))+geom_boxplot(aes(x=gh.species, y=activity.avg))+theme_cowplot()+xlab("Grasshopper species")+ylab("Ant bodyguards (no.)")+theme(legend.position="none")+scale_x_discrete(labels=c("*Orphulella<br>concinnula*", "*Orphulella<br> fluvialis*", "*Orphulella<br>punctata*", "*Orphulella*<br>sp."))+theme(axis.text.x = ggtext::element_markdown())+annotate(geom = "text", x=2.5, y=30, label=paste0("p = ", round(Anova(lm4)[[4]][[5]], 3)))
p10

p11 <- plot_grid(p5, p6, p7, p8, p9, p10, nrow=3, ncol=2, align = "h", labels="AUTO")
p11

save_plot("FigureS1.pdf", p11, base_width=8, base_height=12)
```
