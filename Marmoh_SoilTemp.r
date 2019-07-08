#Set working directory
setwd("H:/US Conservation Program/Germination_Experiments/Marshallia_mohrii&Solidago_ouachitensis/Marshallia mohrii - Germination Data Analysis/2018")

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(lubridate)
library(cowplot)

#Read data
dat = read.csv("Marmoh_SoilTemp.csv",header=T)

#Separate Date.Time column into Date and Time columns
dat2 = dat %>% separate(Date.Time, c("Date", "Time"), sep = " ")
dat2$Date2 = as.Date(dat2$Date,"%m/%d/%Y")
dat2$Month = months(dat2$Date2)
dat2$Year = year(dat2$Date2)

#Subset of only soil data
dat.soil = dat2[dat2$Medium == "Soil",]

#group by hour and day so that temperatures for each time point are averaged between "downslope" and "upslope" sensors
hourly = group_by(dat.soil,Month,Date2,Time) %>%
summarise(HMean = mean(Temp))

#group by day using previous summary to get daily mean, maximum, and minimum from the averages of the two sensors
daily = group_by(hourly,Month,Date2) %>%
summarise(DMean = mean(HMean),
DHigh = max(HMean),
DLow = min(HMean))

#Now get three columns for Min, Max, Mean into one column (longest form possible)
daily.2 = melt(daily, id = c("Month","Date2"))


daily.2$Month <- factor(daily.2$Month, levels = c("March","April","May","June","July","August","September","October","November","December","January","February"))
levels(daily.2$variable)[1] <- 'Mean'
levels(daily.2$variable)[2] <- 'High'
levels(daily.2$variable)[3] <- 'Low'
breaks = c("3/25/2013","5/25/2013","7/25/2013","9/25/2013","11/25/2013","1/25/2014","3/25/2014")
breaks = as.Date(breaks,"%m/%d/%Y")
dbreak = as.Date(daily.2$Date2)[which(unique(as.Date(daily.2$Date2)) %in% breaks)]

#plot
ga = ggplot(data = daily.2,aes(x = Month, y = value, group = variable)) + theme_classic()+
geom_line(aes(x=as.Date(Date2),linetype = variable),lwd = 0.9)+
scale_x_date(breaks = as.Date(dbreak))+
labs(y=expression("Soil Temperature " (degree~C)), x = "Date")+
guides(linetype=guide_legend(title=""))+
theme(panel.grid = element_blank(), panel.border = element_blank(), text = element_text(size = 18), legend.text=element_text(size=14),legend.title=element_text(size=16),
	axis.line=element_line(size=1),axis.text.y=element_text(size=16),axis.text.x=element_text(angle=50,hjust=1,size=16),axis.title=element_text(size=18),
	plot.title = element_text(size=16), legend.position = "top")
ggsave("Marmoh_SoilTemp_Day.tiff",units="in",width=8,height=6,dpi=300)

#Now for month summarize to get the monthly mean temp, the mean high, and the mean minimum temps by averaging across days. 
monthly = group_by(daily,Month) %>%
summarise(Mean = mean(DMean),
SE.mean = sd(Mean)/sqrt(length(Mean)),
High = mean(DHigh),
Low = mean(DLow))

monthly.2 = melt(monthly, id = "Month")

monthly.2$Month <- factor(monthly.2$Month, levels = c("March","April","May","June","July","August","September","October","November","December","January","February"))

gb = ggplot(data = monthly.2[monthly.2$variable != "SE.mean",],aes(x = Month, y = value, group = variable)) + theme_classic()+
geom_line(aes(linetype = variable),lwd = 1)+
labs(y=expression("Soil Temperature " (degree~C)), x = "Month")+
guides(linetype=guide_legend(title=""))+
theme(panel.grid = element_blank(), panel.border = element_blank(), text = element_text(size = 18), legend.text=element_text(size=14),legend.title=element_text(size=16),
	axis.line=element_line(size=1),axis.text.y=element_text(size=16),axis.text.x=element_text(angle=50,hjust=1,size=16),axis.title=element_text(size=18),
	plot.title = element_text(size=16), legend.position = "top")+
ggsave("Marmoh_SoilTemp_Month.tiff",units="in",width=8,height=6,dpi=300)


plot_grid(ga,gb,labels=c("(A)","(B)"),nrow=2)
ggsave("Marmoh_SoilTemp_Comb.tiff",units="in",width=10,height=8,dpi=300)
