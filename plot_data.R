#simple script for plotting data
library(foreign)

#Daphnia cannot live longer than 200 days (in our data not even longer than 123)
lifetime_threshold = 150

#Initialize the variables
age_axis <- c(0:lifetime_threshold)

#Reading the data
dataset <- read.spss("/home/petr/Documents/Master_thesis/AgeEffectsVirulence.sav", to.data.frame = TRUE, use.value.labels = FALSE)
parasite <- as.character(dataset$Pasteuria)
infecteds <- as.numeric(dataset$Infected)

#subsetting
control_population <- dataset[parasite == "Control ",] # Careful, the value of the variable contains the space after Control, possibly FIX the dataset later
host_longevity_control <- control_population$HostLongevity

infected_population <- dataset[infecteds == 1,]
host_longevity_infecteds <-infected_population$HostLongevity

exposed_population <- dataset[(infecteds == 0) & (parasite != "Control "),]
host_longevity_exposed <- exposed_population$HostLongevity


#Survivals of the control, infected and exposed populations
survivals_control <- c(rep(0,lifetime_threshold+1))
survivals_infecteds <- c(rep(0,lifetime_threshold+1))
survivals_exposed <- c(rep(0,lifetime_threshold+1))

for (x in 0:lifetime_threshold) {
  survivals_control[x+1] <- sum(host_longevity_control>x)/ length(host_longevity_control)
  
}

for (x in 0:lifetime_threshold) {
  survivals_infecteds[x+1] <- sum(host_longevity_infecteds>x)/ length(host_longevity_infecteds)
  
}

for (x in 0:lifetime_threshold) {
  survivals_exposed[x+1] <- sum(host_longevity_exposed>x)/ length(host_longevity_exposed)
  
}

pdf("Control_population_bare_survival.pdf")
plot_control_survival <- plot(age_axis, survivals_control, col="blue", xlab = "Age [days]", ylab = "Survivals in control population [fraction]")
legend(0,0.3, legend = c("Data"), col = c("blue"), lty = 1, cex = 0.8)
dev.off()

pdf("Infected_population_bare_survival.pdf")
plot_infecteds_survival <- plot(age_axis, survivals_infecteds, col="blue", xlab = "Age [days]", ylab = "Survivals in infected population [fraction]")
legend(0,0.3, legend = c("Data"), col = c("blue"), lty = 1, cex = 0.8)
dev.off()

pdf("Exposed_population_bare_survival.pdf")
plot_exposed_survival <- plot(age_axis, survivals_exposed, col="blue", xlab = "Age [days]", ylab = "Survivals in exposed population [fraction]")
legend(0,0.3, legend = c("Data"), col = c("blue"), lty = 1, cex = 0.8)
dev.off()

rm(list = ls())