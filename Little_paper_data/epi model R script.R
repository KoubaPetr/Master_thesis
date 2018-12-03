#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------

#OFFSPRING TRAIT MEASUREMENTS

require(cowplot)
#cowplot allows the stictching off ggplot plots together. 

#body size ~ clutch 

read.csv("main.body size.csv")
main.body.size<-read.csv("main.body size.csv")
main.body.size<-main.body.size[-c(36,37,38,39),]#taking out clutch 10 from the data set
main.body.size$Clutch<-as.factor(main.body.size$Clutch)
str(main.body.size$Clutch)
require(ggplot2)
#dot plot
main.bs.se<-as.data.frame(as.list(aggregate(body.size ~ Clutch, data = main.body.size, FUN=function(x) c(mean = mean(x), standerror = sd(x)/sqrt(length(x))))))
#FUN = is how to tell R that you are expecting it to do something to the formula you have written. So 
#in this FUN=function(x) is saying do both following calculations and produce them in two new columns. 

main.bs.plot<-ggplot(data=main.body.size,aes(x=Clutch,y=body.size))+geom_point(position = position_dodge(w=0.75),stat="summary", fun.y="mean")+
  geom_errorbar(data=main.bs.se, aes(x=Clutch, ymin=body.size.mean-1.96*body.size.standerror,ymax=body.size.mean+1.96*body.size.standerror), inherit.aes = FALSE, position=position_dodge(w=0.75))
main.bs.plot<-main.bs.plot + xlab("Clutch") + ylab("Body Size (mm)")
main.bs.plot

#--------------------------------------------------------------------------------------------------------------------------------------

#total reproduction 

read.csv("main.reproduction.csv")
main.reproduction<-read.csv("main.reproduction.csv")
main.reproduction$Clutch<-as.factor(main.reproduction$Clutch)
str(main.reproduction$Clutch)
#dot plot
main.repro.se<-as.data.frame(as.list(aggregate(total.repro ~ Clutch, data = main.reproduction, FUN=function(x) c(mean = mean(x), standerror = sd(x)/sqrt(length(x))))))
#FUN = is how to tell R that you are expecting it to do something to the formula you have written. So 
#in this FUN=function(x) is saying do both following calculations and produce them in two new columns. 

main.repro.plot<-ggplot(data=main.reproduction,aes(x=Clutch,y=total.repro))+geom_point(position = position_dodge(w=0.75),stat="summary", fun.y="mean")+
  geom_errorbar(data=main.repro.se, aes(x=Clutch, ymin=total.repro.mean-1.96*total.repro.standerror,ymax=total.repro.mean+1.96*total.repro.standerror), inherit.aes = FALSE, position=position_dodge(w=0.75))
main.repro.plot<-main.repro.plot + xlab("Clutch") + ylab("Total Offspring")
main.repro.plot

#--------------------------------------------------------------------------------------------------------------------------------------

#infection 

read.csv("main.exposed.csv")
main.exposed<-read.csv("main.exposed.csv")
main.exposed$Clutch<-as.factor(main.exposed$Clutch)
main.exposed<-main.exposed[-c(36,37,38,39),]
str(main.exposed$Clutch)

legend.title<-("Infection Status")#in order to set the legend title, make it an object, then put it in to the manual scale fill. 
main.infect.plot<- ggplot(main.exposed,aes(x = Clutch, fill = I)) + 
  geom_bar(position = "fill") + scale_fill_manual(values=c("black","grey"),legend.title)
main.infect.plot<-main.infect.plot + xlab("Clutch") + ylab("Proportion Infected")
main.infect.plot


#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------

#supplementary data "independently run replication of experiment"
#data in jeff folder
require(cowplot)
#cowplot allows the stictching off ggplot plots together. 

#body size ~ clutch 

read.csv("sm.body size.csv")
sm.body.size<-read.csv("sm.body size.csv")
require(ggplot2)
sm.body.size$Clutch<- as.factor(sm.body.size$Clutch)
#dot plot
sm.bs.se<-as.data.frame(as.list(aggregate(Length ~ Clutch, data = sm.body.size, FUN=function(x) c(mean = mean(x), standerror = sd(x)/sqrt(length(x))))))
#FUN = is how to tell R that you are expecting it to do something to the formula you have written. So 
#in this FUN=function(x) is saying do both following calculations and produce them in two new columns. 

sm.bs.plot<-ggplot(data=sm.body.size,aes(x=Clutch,y=Length))+geom_point(position = position_dodge(w=0.75),stat="summary", fun.y="mean")+
  geom_errorbar(data=sm.bs.se, aes(x=Clutch, ymin=Length.mean-1.96*Length.standerror,ymax=Length.mean+1.96*Length.standerror), inherit.aes = FALSE, position=position_dodge(w=0.75))
sm.bs.plot<-sm.bs.plot + xlab("Clutch") + ylab("Body Size (mm)")
sm.bs.plot

#--------------------------------------------------------------------------------------------------------------------------------------

#total reproduction 

sm.read.csv("sm.total babies.csv")
sm.reproduction<-read.csv("sm.total babies.csv")
sm.reproduction$Clutch<-as.factor(sm.reproduction$Clutch)
#dot plot
sm.repro.se<-as.data.frame(as.list(aggregate(Total.Baby.Number ~ Clutch, data = sm.reproduction, FUN=function(x) c(mean = mean(x), standerror = sd(x)/sqrt(length(x))))))
#FUN = is how to tell R that you are expecting it to do something to the formula you have written. So 
#in this FUN=function(x) is saying do both following calculations and produce them in two new columns. 

sm.repro.plot<-ggplot(data=sm.reproduction,aes(x=Clutch,y=Total.Baby.Number))+geom_point(position = position_dodge(w=0.75),stat="summary", fun.y="mean")+
  geom_errorbar(data=sm.repro.se, aes(x=Clutch, ymin=Total.Baby.Number.mean-1.96*Total.Baby.Number.standerror,ymax=Total.Baby.Number.mean+1.96*Total.Baby.Number.standerror), inherit.aes = FALSE, position=position_dodge(w=0.75))
sm.repro.plot<-sm.repro.plot + xlab("Clutch") + ylab("Total Offspring")
sm.repro.plot

#--------------------------------------------------------------------------------------------------------------------------------------

#infection 

  read.csv("sm.infection status.csv")
  sm.exposed<-read.csv("sm.infection status.csv")
  sm.exposed<-subset(sm.exposed, Exposed. == "e")
  sm.exposed$Clutch<-as.factor(sm.exposed$Clutch)
  str(sm.exposed$Clutch)
  legend.title<-"Infection Status"#in order to set the legend title, make it an object, then put it in to the manual scale fill. 
  sm.infect.plot<- ggplot(sm.exposed,aes(x = Clutch, fill = Infected.)) + 
    geom_bar(position = "fill") + scale_fill_manual(values=c("black","grey"),legend.title)
  sm.infect.plot<-sm.infect.plot + xlab("Clutch") + ylab("Proportion Infected")
  sm.infect.plot
#--------------------------------------------------------------------------------------------------------------------------------------

#combine plots into a panel

epi.model.add.plots<-plot_grid(bs.plot, infect.plot, repro.plot, labels=c("A", "B", "C"), ncol = 1, nrow = 3)
save_plot("epi.model.add.plots.pdf", epi.model.add.plots, base_aspect_ratio = 1.3)

#------------------------------------------------------------------------------------------------------------------------

#suggested put both the sm and main plots over one another into one panel. 

#body size

all.body.size<- ggplot()+
  geom_point(data=main.body.size,aes(x=Clutch,y=body.size),position = position_dodge(w=0.75),stat="summary", fun.y="mean")+
  geom_errorbar(data=main.bs.se, aes(x=Clutch, ymin=body.size.mean-1.96*body.size.standerror,ymax=body.size.mean+1.96*body.size.standerror),width = 0.1, inherit.aes = FALSE, position=position_dodge(w=0.75)) +
  geom_point(data=sm.body.size, aes(x=Clutch,y=Length), position = position_dodge(w=0.75),stat="summary", fun.y="mean", color="grey")+
  geom_errorbar(data=sm.bs.se, aes(x=Clutch, ymin=Length.mean-1.96*Length.standerror,ymax=Length.mean+1.96*Length.standerror), width = 0.1, color = "grey", inherit.aes = FALSE, position=position_dodge(w=0.75))+

  xlab("Clutch") + ylab("Offspring Body Size (mm)")
all.body.size

#offspring/ total reproduction
all.repro<- ggplot()+
  geom_point(data=main.reproduction,aes(x=Clutch,y=total.repro), position = position_dodge(w=0.75),stat="summary", fun.y="mean")+
  geom_errorbar(data=main.repro.se, aes(x=Clutch, ymin=total.repro.mean-1.96*total.repro.standerror,ymax=total.repro.mean+1.96*total.repro.standerror), width = 0.1, inherit.aes = FALSE, position=position_dodge(w=0.75)) +
  geom_point(data=sm.reproduction, aes(x=Clutch,y=Total.Baby.Number),position = position_dodge(w=0.75),stat="summary", fun.y="mean", color="grey") +
  geom_errorbar(data=sm.repro.se, aes(x=Clutch, ymin=Total.Baby.Number.mean-1.96*Total.Baby.Number.standerror,ymax=Total.Baby.Number.mean+1.96*Total.Baby.Number.standerror),width = 0.1, color = "grey", inherit.aes = FALSE, position=position_dodge(w=0.75))+ 

  xlab("Clutch") + ylab("Total Reproduction")
all.repro

#-------

#infection

combine.infect<-read.csv("combined.infect.proportions.csv")
str(combine.infect)
combine.infect$clutch<-as.factor(combine.infect$clutch)
combine.infect$experiment<-as.factor(combine.infect$experiment)
combine.infect.plot<-ggplot(combine.infect)+
  geom_bar(aes(x=clutch, y=prop.infect, fill=experiment, order=clutch) ,
           stat="identity",
           position = "dodge")+ 
  geom_errorbar(data=combine.infect, aes(x=clutch, ymin=lower.confidence.levels,ymax=upper.confidence.levels),width = 0.1, inherit.aes = FALSE, position=position_dodge(w=0.77))+
  ylab("proportion infected")+ scale_fill_manual(values=c("black","grey"))
combine.infect.plot
#------
combine.infect.exp1<-subset(combine.infect, experiment == 1)

combine.infect.exp2<-subset(combine.infect, experiment == 2)

combine.infect.plot<- ggplot()+
  geom_point(data=combine.infect.exp1,aes(x=clutch,y=prop.infect), position = position_dodge(w=0.75),stat="summary", fun.y="mean")+
  geom_errorbar(data=combine.infect.exp1, aes(x=clutch, ymin=lower.confidence.levels,ymax=upper.confidence.levels), width = 0.1, inherit.aes = FALSE, position=position_dodge(w=0.75)) +
  geom_point(data=combine.infect.exp2, aes(x=clutch,y=prop.infect),position = position_dodge(w=0.75),stat="summary", fun.y="mean", color="grey") +
 geom_errorbar(data=combine.infect.exp2, aes(x=clutch, ymin=lower.confidence.levels,ymax=upper.confidence.levels),width = 0.1, color = "grey", inherit.aes = FALSE, position=position_dodge(w=0.75))+
  xlab("Clutch") + ylab("Proportion Infected")
combine.infect.plot



require(cowplot)
all.exp.plots.combine<-plot_grid(combine.infect.plot, all.repro, all.body.size, labels=c("A", "B", "C"), ncol = 1, nrow = 3)

all.exp.plots.combine
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
#NO PATHOGEN PRESENT

# all rates are per day
r <- 3 # max clutch size/number of days between clutches
K <- 1 # max pop density
a <- 1/5 # additional parasite induced mortality
m <- 1/10 # maturation rate (average age at maturation is 10 days)
BYy <- 3.5 # transmission to young individuals with young mothers
mat1<-0 # maternal effect on susceptibility
age1<-0 # ageing effect on susceptibility
BYo <- BYy*(1-mat1) # transmission to young individuals with old mothers
BOy <- BYy*(1-age1) # transmission to old individuals with young mothers
BOo <- BYy*(1-mat1)*(1-age1) # transmission to old individuals with old mothers
p<-0 # proportion of natural deaths leading to transmission


ds<-seq(from=1/1000, to=1/1,length.out=101) # this is the list of values of death rate
no_path_eq<-matrix(nrow=length(ds),ncol=4) # to store equilibrium with pathogen absent
w_path_eq<-matrix(nrow=length(ds),ncol=8) # to store equilibrium with pathogen


# equilibrium with no pathogen (I've analytically worked out the equilibria using Mathematica)
no_path_eq[,1]<-((ds^2)*K*r-(ds^3)*K)/(r*((ds+m)^2)) #UYy
no_path_eq[,2]<-((r-ds)*ds*K*m)/(r*((ds+m)^2))  #UYo
no_path_eq[,3]<-((r-ds)*ds*K*m)/(r*((ds+m)^2)) #UOy
no_path_eq[,4]<-((r-ds)*K*(m^2))/(r*((ds+m)^2)) #UOo

#this is the plotting of the population dynamics with no pathogen present. 
 
dev.new()
cols<-c("firebrick1","dodgerblue1","goldenrod1","forestgreen", "blueviolet") # colours for plots
plot(NA,NA,xlim=c(0,max(ds)),ylim=c(0,1),xlab="Mortality", ylab = "Population Density")

lines(no_path_eq[,1]~ds,lwd=2, lty = 1,col="firebrick1")

lines(no_path_eq[,2]~ds,lwd=5, lty = 1,col="dodgerblue1")

lines(no_path_eq[,3]~ds,lwd=2, lty=1,col="goldenrod1")

lines(no_path_eq[,4]~ds,lwd=2, lty=1,col="forestgreen")

lines(no_path_eq[,1]+no_path_eq[,2]+no_path_eq[,3]+no_path_eq[,4]~ds, lwd=2, lty=1,col="blueviolet")

#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------


#model output


library(deSolve)
library(rootSolve)

# Yyu refers to young individual who had a young mother and is uninfected...
# Oyi refers to an old individual who had a young mother and is infected... etc.
# r is the growth rate (maximum clutch size)
# K is the nutrient limitation carrying capacity (population density at which clutch size goes to zero)
# d is the baseline mortality/death rate
# m is the maturation rate (rate at which individuals switch from young to old)
# Bxx is the transmission rate to an uninfected individual of type xx, BYy = yound individual with young mother, BOy = old individual with young mother, etc.
# assume infected individuals don't reproduce
# initially ignore aging (assume d is the same for young and old)


#PLOTS P = 0
#plot NO effects
# all rates are per day
r <- 3 # max clutch size/number of days between clutches
K <- 1 # max pop density
#a <- 1/10 # additional parasite induced mortality
#m<- 1/10 # maturation rate (average age at maturation is 10 days)
BYy <- 3.5 # transmission to young individuals with young mothers
mat1<-0 # maternal effect on susceptibility
age1<-0 # ageing effect on susceptibility
BYo <- BYy*(1-mat1) # transmission to young individuals with old mothers
BOy <- BYy*(1-age1) # transmission to old individuals with young mothers
BOo <- BYy*(1-mat1)*(1-age1) # transmission to old individuals with old mothers
p<-0 # proportion of natural deaths leading to transmission
#set to zero would be no deaths leading to transmission because there would be no disease

ds<-seq(from=1/1000, to=1/1,length.out=101) # this is the list of values of death rate
no_path_eq<-matrix(nrow=length(ds),ncol=4) # to store equilibrium with pathogen absent
w_path_eq<-matrix(nrow=length(ds),ncol=8) # to store equilibrium with pathogen

# equilibrium with no pathogen (I've analytically worked out the equilibria using Mathematica)
no_path_eq[,1]<-((ds^2)*K*r-(ds^3)*K)/(r*((ds+m)^2)) #UYy
no_path_eq[,2]<-((r-ds)*ds*K*m)/(r*((ds+m)^2))  #UYo
no_path_eq[,3]<-((r-ds)*ds*K*m)/(r*((ds+m)^2)) #UOy
no_path_eq[,4]<-((r-ds)*K*(m^2))/(r*((ds+m)^2)) #UOo



weights<-matrix(ncol=4,nrow=length(ds)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums
trans_pot1<-rowSums(no_path_eq*weights)*(a+ds*p)/(a+ds) # caculate transmission potential. This is sum of the abundance of individuals of each class multiplied by their susceptibility
#-----------------------------------------------------------------------------

#plot ONE effect
# all rates are per day
r <- 3 # max clutch size/number of days between clutches
K <- 1 # max pop density
#a <- 1/10 # additional parasite induced mortality
#m<- 1/10 # maturation rate (average age at maturation is 10 days)
BYy <- 3.5 # transmission to young individuals with young mothers
mat2<-0.5 # maternal effect on susceptibility
age2<-0 # ageing effect on susceptibility
BYo <- BYy*(1-mat2) # transmission to young individuals with old mothers
BOy <- BYy*(1-age2) # transmission to old individuals with young mothers
BOo <- BYy*(1-mat2)*(1-age2) # transmission to old individuals with old mothers
p<-0 # proportion of natural deaths leading to transmission
#set to zero would be no deaths leading to transmission because there would be no disease

#ds<-seq(from=1/1000, to=1/5,length.out=101) # this is the list of values of death rate
no_path_eq<-matrix(nrow=length(ds),ncol=4) # to store equilibrium with pathogen absent
w_path_eq<-matrix(nrow=length(ds),ncol=8) # to store equilibrium with pathogen

# equilibrium with no pathogen (I've analytically worked out the equilibria using Mathematica)
no_path_eq[,1]<-((ds^2)*K*r-(ds^3)*K)/(r*((ds+m)^2)) #UYy
no_path_eq[,2]<-((r-ds)*ds*K*m)/(r*((ds+m)^2))  #UYo
no_path_eq[,3]<-((r-ds)*ds*K*m)/(r*((ds+m)^2)) #UOy
no_path_eq[,4]<-((r-ds)*K*(m^2))/(r*((ds+m)^2)) #UOo



weights<-matrix(ncol=4,nrow=length(ds)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums
trans_pot2<-rowSums(no_path_eq*weights)*(a+ds*p)/(a+ds) # caculate transmission potential. This is sum of the abundance of individuals of each class multiplied by their susceptibility

#plot ONE effect
# all rates are per day
r <- 3 # max clutch size/number of days between clutches
K <- 1 # max pop density
#a <- 1/10 # additional parasite induced mortality
#m<- 1/10 # maturation rate (average age at maturation is 10 days)
BYy <- 3.5 # transmission to young individuals with young mothers
mat2a<-0 # maternal effect on susceptibility
age2a<-0.5 # ageing effect on susceptibility
BYo <- BYy*(1-mat2a) # transmission to young individuals with old mothers
BOy <- BYy*(1-age2a) # transmission to old individuals with young mothers
BOo <- BYy*(1-mat2a)*(1-age2a) # transmission to old individuals with old mothers
p<-0 # proportion of natural deaths leading to transmission
#set to zero would be no deaths leading to transmission because there would be no disease

#ds<-seq(from=1/1000, to=1/5,length.out=101) # this is the list of values of death rate
no_path_eq<-matrix(nrow=length(ds),ncol=4) # to store equilibrium with pathogen absent
w_path_eq<-matrix(nrow=length(ds),ncol=8) # to store equilibrium with pathogen

# equilibrium with no pathogen (I've analytically worked out the equilibria using Mathematica)
no_path_eq[,1]<-((ds^2)*K*r-(ds^3)*K)/(r*((ds+m)^2)) #UYy
no_path_eq[,2]<-((r-ds)*ds*K*m)/(r*((ds+m)^2))  #UYo
no_path_eq[,3]<-((r-ds)*ds*K*m)/(r*((ds+m)^2)) #UOy
no_path_eq[,4]<-((r-ds)*K*(m^2))/(r*((ds+m)^2)) #UOo



weights<-matrix(ncol=4,nrow=length(ds)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums
trans_pot2A<-rowSums(no_path_eq*weights)*(a+ds*p)/(a+ds) # caculate transmission potential. This is sum of the abundance of individuals of each class multiplied by their susceptibility
#----------
#plot ALL effects
# all rates are per day
r <- 3 # max clutch size/number of days between clutches
K <- 1 # max pop density
#a <- 1/10 # additional parasite induced mortality
#m<- 1/10 # maturation rate (average age at maturation is 10 days)
BYy <- 3.5 # transmission to young individuals with young mothers
mat3<-0.5 # maternal effect on susceptibility
age3<-0.5 # ageing effect on susceptibility
BYo <- BYy*(1-mat3) # transmission to young individuals with old mothers
BOy <- BYy*(1-age3) # transmission to old individuals with young mothers
BOo <- BYy*(1-mat3)*(1-age3) # transmission to old individuals with old mothers
p<-0 # proportion of natural deaths leading to transmission
#set to zero would be no deaths leading to transmission because there would be no disease

#ds<-seq(from=1/1000, to=1/5,length.out=101) # this is the list of values of death rate
no_path_eq<-matrix(nrow=length(ds),ncol=4) # to store equilibrium with pathogen absent
w_path_eq<-matrix(nrow=length(ds),ncol=8) # to store equilibrium with pathogen

# equilibrium with no pathogen (I've analytically worked out the equilibria using Mathematica)
no_path_eq[,1]<-((ds^2)*K*r-(ds^3)*K)/(r*((ds+m)^2)) #UYy
no_path_eq[,2]<-((r-ds)*ds*K*m)/(r*((ds+m)^2))  #UYo
no_path_eq[,3]<-((r-ds)*ds*K*m)/(r*((ds+m)^2)) #UOy
no_path_eq[,4]<-((r-ds)*K*(m^2))/(r*((ds+m)^2)) #UOo



weights<-matrix(ncol=4,nrow=length(ds)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums
trans_pot3<-rowSums(no_path_eq*weights)*(a+ds*p)/(a+ds) # caculate transmission potential. This is sum of the abundance of individuals of each class multiplied by their susceptibility
#----------
#plot P = 0.5 effects
# all rates are per day
r <- 3 # max clutch size/number of days between clutches
K <- 1 # max pop density
#a <- 1/10 # additional parasite induced mortality
#m<- 1/10 # maturation rate (average age at maturation is 10 days)
BYy <- 3.5 # transmission to young individuals with young mothers
mat4<-0 # maternal effect on susceptibility
age4<-0 # ageing effect on susceptibility
BYo <- BYy*(1-mat4) # transmission to young individuals with old mothers
BOy <- BYy*(1-age4) # transmission to old individuals with young mothers
BOo <- BYy*(1-mat4)*(1-age4) # transmission to old individuals with old mothers
p<-0.5 # proportion of natural deaths leading to transmission
#set to zero would be no deaths leading to transmission because there would be no disease

#ds<-seq(from=1/1000, to=1/5,length.out=101) # this is the list of values of death rate
no_path_eq<-matrix(nrow=length(ds),ncol=4) # to store equilibrium with pathogen absent
w_path_eq<-matrix(nrow=length(ds),ncol=8) # to store equilibrium with pathogen

# equilibrium with no pathogen (I've analytically worked out the equilibria using Mathematica)
no_path_eq[,1]<-((ds^2)*K*r-(ds^3)*K)/(r*((ds+m)^2)) #UYy
no_path_eq[,2]<-((r-ds)*ds*K*m)/(r*((ds+m)^2))  #UYo
no_path_eq[,3]<-((r-ds)*ds*K*m)/(r*((ds+m)^2)) #UOy
no_path_eq[,4]<-((r-ds)*K*(m^2))/(r*((ds+m)^2)) #UOo



weights<-matrix(ncol=4,nrow=length(ds)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums
trans_pot4<-rowSums(no_path_eq*weights)*(a+ds*p)/(a+ds) # caculate transmission potential. This is sum of the abundance of individuals of each class multiplied by their susceptibility
#-----------------------------------------------------------------------------

#plot ONE effect
# all rates are per day
r <- 3 # max clutch size/number of days between clutches
K <- 1 # max pop density
#a <- 1/10 # additional parasite induced mortality
#m<- 1/10 # maturation rate (average age at maturation is 10 days)
BYy <- 3.5 # transmission to young individuals with young mothers
mat5<-0.5 # maternal effect on susceptibility
age5<-0 # ageing effect on susceptibility
BYo <- BYy*(1-mat5) # transmission to young individuals with old mothers
BOy <- BYy*(1-age5) # transmission to old individuals with young mothers
BOo <- BYy*(1-mat5)*(1-age5) # transmission to old individuals with old mothers
p<-0.5 # proportion of natural deaths leading to transmission
#set to zero would be no deaths leading to transmission because there would be no disease

#ds<-seq(from=1/1000, to=1/5,length.out=101) # this is the list of values of death rate
no_path_eq<-matrix(nrow=length(ds),ncol=4) # to store equilibrium with pathogen absent
w_path_eq<-matrix(nrow=length(ds),ncol=8) # to store equilibrium with pathogen

# equilibrium with no pathogen (I've analytically worked out the equilibria using Mathematica)
no_path_eq[,1]<-((ds^2)*K*r-(ds^3)*K)/(r*((ds+m)^2)) #UYy
no_path_eq[,2]<-((r-ds)*ds*K*m)/(r*((ds+m)^2))  #UYo
no_path_eq[,3]<-((r-ds)*ds*K*m)/(r*((ds+m)^2)) #UOy
no_path_eq[,4]<-((r-ds)*K*(m^2))/(r*((ds+m)^2)) #UOo



weights<-matrix(ncol=4,nrow=length(ds)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums
trans_pot5<-rowSums(no_path_eq*weights)*(a+ds*p)/(a+ds) # caculate transmission potential. This is sum of the abundance of individuals of each class multiplied by their susceptibility

#plot ONE effect
# all rates are per day
r <- 3 # max clutch size/number of days between clutches
K <- 1 # max pop density
#a <- 1/10 # additional parasite induced mortality
#m<- 1/10 # maturation rate (average age at maturation is 10 days)
BYy <- 3.5 # transmission to young individuals with young mothers
mat5a<-0 # maternal effect on susceptibility
age5a<-0.5 # ageing effect on susceptibility
BYo <- BYy*(1-mat5a) # transmission to young individuals with old mothers
BOy <- BYy*(1-age5a) # transmission to old individuals with young mothers
BOo <- BYy*(1-mat5a)*(1-age5a) # transmission to old individuals with old mothers
p<-0.5 # proportion of natural deaths leading to transmission
#set to zero would be no deaths leading to transmission because there would be no disease

#ds<-seq(from=1/1000, to=1/5,length.out=101) # this is the list of values of death rate
no_path_eq<-matrix(nrow=length(ds),ncol=4) # to store equilibrium with pathogen absent
w_path_eq<-matrix(nrow=length(ds),ncol=8) # to store equilibrium with pathogen

# equilibrium with no pathogen (I've analytically worked out the equilibria using Mathematica)
no_path_eq[,1]<-((ds^2)*K*r-(ds^3)*K)/(r*((ds+m)^2)) #UYy
no_path_eq[,2]<-((r-ds)*ds*K*m)/(r*((ds+m)^2))  #UYo
no_path_eq[,3]<-((r-ds)*ds*K*m)/(r*((ds+m)^2)) #UOy
no_path_eq[,4]<-((r-ds)*K*(m^2))/(r*((ds+m)^2)) #UOo



weights<-matrix(ncol=4,nrow=length(ds)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums
trans_pot5a<-rowSums(no_path_eq*weights)*(a+ds*p)/(a+ds) # caculate transmission potential. This is sum of the abundance of individuals of each class multiplied by their susceptibility
#----------

#plot ALL effects
# all rates are per day
r <- 3 # max clutch size/number of days between clutches
K <- 1 # max pop density
#a <- 1/10 # additional parasite induced mortality
#m<- 1/10 # maturation rate (average age at maturation is 10 days)
BYy <- 3.5 # transmission to young individuals with young mothers
mat6<-0.5 # maternal effect on susceptibility
age6<-0.5 # ageing effect on susceptibility
BYo <- BYy*(1-mat6) # transmission to young individuals with old mothers
BOy <- BYy*(1-age6) # transmission to old individuals with young mothers
BOo <- BYy*(1-mat6)*(1-age6) # transmission to old individuals with old mothers
p<-0.5 # proportion of natural deaths leading to transmission
#set to zero would be no deaths leading to transmission because there would be no disease

#ds<-seq(from=1/1000, to=1/5,length.out=101) # this is the list of values of death rate
no_path_eq<-matrix(nrow=length(ds),ncol=4) # to store equilibrium with pathogen absent
w_path_eq<-matrix(nrow=length(ds),ncol=8) # to store equilibrium with pathogen

# equilibrium with no pathogen (I've analytically worked out the equilibria using Mathematica)
no_path_eq[,1]<-((ds^2)*K*r-(ds^3)*K)/(r*((ds+m)^2)) #UYy
no_path_eq[,2]<-((r-ds)*ds*K*m)/(r*((ds+m)^2))  #UYo
no_path_eq[,3]<-((r-ds)*ds*K*m)/(r*((ds+m)^2)) #UOy
no_path_eq[,4]<-((r-ds)*K*(m^2))/(r*((ds+m)^2)) #UOo



weights<-matrix(ncol=4,nrow=length(ds)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums
trans_pot6<-rowSums(no_path_eq*weights)*(a+ds*p)/(a+ds) # caculate transmission potential. This is sum of the abundance of individuals of each class multiplied by their susceptibility
#----------
#P = 1
#plot NO effects
# all rates are per day
r <- 3 # max clutch size/number of days between clutches
K <- 1 # max pop density
#a <- 1/10 # additional parasite induced mortality
#m<- 1/10 # maturation rate (average age at maturation is 10 days)
BYy <- 3.5 # transmission to young individuals with young mothers
mat7<-0 # maternal effect on susceptibility
age7<-0 # ageing effect on susceptibility
BYo <- BYy*(1-mat7) # transmission to young individuals with old mothers
BOy <- BYy*(1-age7) # transmission to old individuals with young mothers
BOo <- BYy*(1-mat7)*(1-age7) # transmission to old individuals with old mothers
p<-1 # proportion of natural deaths leading to transmission
#set to zero would be no deaths leading to transmission because there would be no disease

#ds<-seq(from=1/1000, to=1/5,length.out=101) # this is the list of values of death rate
no_path_eq<-matrix(nrow=length(ds),ncol=4) # to store equilibrium with pathogen absent
w_path_eq<-matrix(nrow=length(ds),ncol=8) # to store equilibrium with pathogen

# equilibrium with no pathogen (I've analytically worked out the equilibria using Mathematica)
no_path_eq[,1]<-((ds^2)*K*r-(ds^3)*K)/(r*((ds+m)^2)) #UYy
no_path_eq[,2]<-((r-ds)*ds*K*m)/(r*((ds+m)^2))  #UYo
no_path_eq[,3]<-((r-ds)*ds*K*m)/(r*((ds+m)^2)) #UOy
no_path_eq[,4]<-((r-ds)*K*(m^2))/(r*((ds+m)^2)) #UOo



weights<-matrix(ncol=4,nrow=length(ds)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums
trans_pot7<-rowSums(no_path_eq*weights)*(a+ds*p)/(a+ds) # caculate transmission potential. This is sum of the abundance of individuals of each class multiplied by their susceptibility
#-----------------------------------------------------------------------------


#plot ONE effect
# all rates are per day
r <- 3 # max clutch size/number of days between clutches
K <- 1 # max pop density
#a <- 1/10 # additional parasite induced mortality
#m<- 1/10 # maturation rate (average age at maturation is 10 days)
BYy <- 3.5 # transmission to young individuals with young mothers
mat8<-0.5 # maternal effect on susceptibility
age8<-0 # ageing effect on susceptibility
BYo <- BYy*(1-mat8) # transmission to young individuals with old mothers
BOy <- BYy*(1-age8) # transmission to old individuals with young mothers
BOo <- BYy*(1-mat8)*(1-age8) # transmission to old individuals with old mothers
p<-1 # proportion of natural deaths leading to transmission
#set to zero would be no deaths leading to transmission because there would be no disease

#ds<-seq(from=1/1000, to=1/5,length.out=101) # this is the list of values of death rate
no_path_eq<-matrix(nrow=length(ds),ncol=4) # to store equilibrium with pathogen absent
w_path_eq<-matrix(nrow=length(ds),ncol=8) # to store equilibrium with pathogen

# equilibrium with no pathogen (I've analytically worked out the equilibria using Mathematica)
no_path_eq[,1]<-((ds^2)*K*r-(ds^3)*K)/(r*((ds+m)^2)) #UYy
no_path_eq[,2]<-((r-ds)*ds*K*m)/(r*((ds+m)^2))  #UYo
no_path_eq[,3]<-((r-ds)*ds*K*m)/(r*((ds+m)^2)) #UOy
no_path_eq[,4]<-((r-ds)*K*(m^2))/(r*((ds+m)^2)) #UOo



weights<-matrix(ncol=4,nrow=length(ds)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums
trans_pot8<-rowSums(no_path_eq*weights)*(a+ds*p)/(a+ds) # caculate transmission potential. This is sum of the abundance of individuals of each class multiplied by their susceptibility

#plot ONE effect
# all rates are per day
r <- 3 # max clutch size/number of days between clutches
K <- 1 # max pop density
#a <- 1/10 # additional parasite induced mortality
#m<- 1/10 # maturation rate (average age at maturation is 10 days)
BYy <- 3.5 # transmission to young individuals with young mothers
mat8a<-0 # maternal effect on susceptibility
age8a<-0.5 # ageing effect on susceptibility
BYo <- BYy*(1-mat8a) # transmission to young individuals with old mothers
BOy <- BYy*(1-age8a) # transmission to old individuals with young mothers
BOo <- BYy*(1-mat8a)*(1-age8a) # transmission to old individuals with old mothers
p<-1 # proportion of natural deaths leading to transmission
#set to zero would be no deaths leading to transmission because there would be no disease

#ds<-seq(from=1/1000, to=1/5,length.out=101) # this is the list of values of death rate
no_path_eq<-matrix(nrow=length(ds),ncol=4) # to store equilibrium with pathogen absent
w_path_eq<-matrix(nrow=length(ds),ncol=8) # to store equilibrium with pathogen

# equilibrium with no pathogen (I've analytically worked out the equilibria using Mathematica)
no_path_eq[,1]<-((ds^2)*K*r-(ds^3)*K)/(r*((ds+m)^2)) #UYy
no_path_eq[,2]<-((r-ds)*ds*K*m)/(r*((ds+m)^2))  #UYo
no_path_eq[,3]<-((r-ds)*ds*K*m)/(r*((ds+m)^2)) #UOy
no_path_eq[,4]<-((r-ds)*K*(m^2))/(r*((ds+m)^2)) #UOo



weights<-matrix(ncol=4,nrow=length(ds)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums
trans_pot8a<-rowSums(no_path_eq*weights)*(a+ds*p)/(a+ds) # caculate transmission potential. This is sum of the abundance of individuals of each class multiplied by their susceptibility
#----------
#----------
#plot ALL effects
# all rates are per day
r <- 3 # max clutch size/number of days between clutches
K <- 1 # max pop density
#a <- 1/10 # additional parasite induced mortality
#m<- 1/10 # maturation rate (average age at maturation is 10 days)
BYy <- 3.5 # transmission to young individuals with young mothers
mat9<-0.5 # maternal effect on susceptibility
age9<-0.5 # ageing effect on susceptibility
BYo <- BYy*(1-mat9) # transmission to young individuals with old mothers
BOy <- BYy*(1-age9) # transmission to old individuals with young mothers
BOo <- BYy*(1-mat9)*(1-age9) # transmission to old individuals with old mothers
p<-1 # proportion of natural deaths leading to transmission
#set to zero would be no deaths leading to transmission because there would be no disease

#ds<-seq(from=1/1000, to=1/5,length.out=101) # this is the list of values of death rate
no_path_eq<-matrix(nrow=length(ds),ncol=4) # to store equilibrium with pathogen absent
w_path_eq<-matrix(nrow=length(ds),ncol=8) # to store equilibrium with pathogen

# equilibrium with no pathogen (I've analytically worked out the equilibria using Mathematica)
no_path_eq[,1]<-((ds^2)*K*r-(ds^3)*K)/(r*((ds+m)^2)) #UYy
no_path_eq[,2]<-((r-ds)*ds*K*m)/(r*((ds+m)^2))  #UYo
no_path_eq[,3]<-((r-ds)*ds*K*m)/(r*((ds+m)^2)) #UOy
no_path_eq[,4]<-((r-ds)*K*(m^2))/(r*((ds+m)^2)) #UOo



weights<-matrix(ncol=4,nrow=length(ds)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums
trans_pot9<-rowSums(no_path_eq*weights)*(a+ds*p)/(a+ds) # caculate transmission potential. This is sum of the abundance of individuals of each class multiplied by their susceptibility
#----------

#FIRST THREE WILL BE P = 0 
#PLOT 1 NO EFFECTS PLOT 
#PLOT 2 ONE EFFECT 
#PLOT 3 BOTH EFFECTS
cols<-c("firebrick1","dodgerblue1","goldenrod1","forestgreen") # colours for plots

dev.new()
par(mfrow=c(2,2)) #mar = margins, figure out the order
#PLOT 1
plot(NA,NA,xlim=c(0,max(ds)),ylim=c(0,3.5),xlab=" Mortality", ylab="Transmission Potential")
title(outer=TRUE,adj=0.1,main="A",cex=1, col="black",font=1,line=-4)
i<-1
lines(trans_pot1~ds,col=cols[i],lwd=1)
i<-2
lines(trans_pot4~ds,col=cols[i],lwd=1)
i<-3
lines(trans_pot7~ds,col=cols[i],lwd=1)


#PLOT 2
plot(NA,NA,xlim=c(0,max(ds)),ylim=c(0,3.5),xlab=" Mortality", ylab="Transmission Potential")
title(outer=TRUE,adj=0.61,main="B",cex=1.1, col="black",font=1,line=-4)
i<-1
lines(trans_pot2~ds,col=cols[i],lwd=1)
i<-2
lines(trans_pot5~ds,col=cols[i],lwd=1)
i<-3
lines(trans_pot8~ds, col=cols[i],lwd=1)

#PLOT 3
plot(NA,NA,xlim=c(0,max(ds)),ylim=c(0,3.5),xlab=" Mortality", ylab="Transmission Potential")
title(outer=TRUE,adj=0.1,main="C",cex=1.1, col="black",font=1,line=-25)
i<-1
lines(trans_pot2A~ds,col=cols[i],lwd=1)
i<-2
lines(trans_pot5a~ds,col=cols[i],lwd=1)
i<-3
lines(trans_pot8a~ds, col=cols[i],lwd=1)


#PLOT 4
plot(NA,NA,xlim=c(0,max(ds)),ylim=c(0,3.5),xlab=" Mortality", ylab="Transmission Potential")
title(outer=TRUE,adj=0.61,main="D",cex=1.1, col="black",font=1,line=-25)
i<-1
lines(trans_pot3~ds,col=cols[i],lwd=1)
i<-2
lines(trans_pot6~ds,col=cols[i],lwd=1)
i<-3
lines(trans_pot9~ds,col=cols[i],lwd=1)







