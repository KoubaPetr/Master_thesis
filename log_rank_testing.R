# log-rank testing

require('survival')

### Example - lung
data(lung)

survival_object <- Surv(time = lung$time, event = lung$status)
mfit <- survfit(survival_object ~ lung$sex)
plot(mfit, col = c("red", "blue"))

survdiff(survival_object ~ lung$sex)

### Our data

require('foreign')

## reading data
dataset <- read.spss("/home/petr/Documents/Master_thesis/AgeEffectsVirulence.sav", to.data.frame = TRUE, use.value.labels = FALSE)
infected_at_5 <- dataset[dataset$InfectionAge == 5 & dataset$Infected == 1,]
infected_at_15 <- dataset[dataset$InfectionAge == 15 & dataset$Infected == 1,]
infected_at_30 <- dataset[dataset$InfectionAge == 30 & dataset$Infected == 1,]

infected_at_5_or_15 <- dataset[dataset$InfectionAge != 30 & dataset$Infected == 1,]
infected_at_15_or_30 <- dataset[dataset$InfectionAge != 5 & dataset$Infected == 1,]
infected_at_5_or_30 <- dataset[dataset$InfectionAge != 15 & dataset$Infected == 1,]


## preparing survival objects
surv_infected_at_5 <- Surv(time = infected_at_5$HostLongevity, event = rep(1, length(infected_at_5$HostLongevity)))
surv_infected_at_15 <- Surv(time = infected_at_15$HostLongevity, event = rep(1, length(infected_at_15$HostLongevity)))
surv_infected_at_30 <- Surv(time = infected_at_30$HostLongevity, event = rep(1, length(infected_at_30$HostLongevity)))

surv_infected_at_5_or_15 <- Surv(time = infected_at_5_or_15$HostLongevity, event = rep(1, length(infected_at_5_or_15$HostLongevity)))
surv_infected_at_15_or_30 <- Surv(time = infected_at_15_or_30$HostLongevity, event = rep(1, length(infected_at_15_or_30$HostLongevity)))
surv_infected_at_5_or_30 <- Surv(time = infected_at_5_or_30$HostLongevity, event = rep(1, length(infected_at_5_or_30$HostLongevity)))

plot(survfit(surv_infected_at_5_or_15 ~ infected_at_5_or_15$InfectionAge), col = c('red', 'blue'))
lines(0:150,lapply(0:150, function(x){return(sum(infected_at_5_or_15$HostLongevity>x)/ sum(infected_at_5_or_15$HostLongevity>15))}), type = 's')
survdiff(surv_infected_at_5_or_15 ~ infected_at_5_or_15$InfectionAge)

plot(survfit(surv_infected_at_15_or_30 ~ infected_at_15_or_30$InfectionAge), col = c('red', 'blue'))
survdiff(surv_infected_at_15_or_30 ~ infected_at_15_or_30$InfectionAge)

plot(survfit(surv_infected_at_5_or_30 ~ infected_at_5_or_30$InfectionAge), col = c('red', 'blue'))
survdiff(surv_infected_at_5_or_30 ~ infected_at_5_or_30$InfectionAge)

# TODO: Why does the survival estimated by survfit() differ from the actual survival?
#       - Is it because the death times are not unique? (It seems that is a prerequisite for survfit() )
#       - Is there some other method which handles that? Or how could we actually use survfit() in our case?
#       - Also try to run the log-rank test for the age of infection longevities, so it would be comparable
#         Maybe ask about this on stack exchange?
#       Try the method using LRT, as discussed with Roland
#       - think of different Ansatzes to be used, so far most of them implemented in delta_inf(), but ansatz a^2 + const,
#         could also be of an interest

