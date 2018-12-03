#Estimation of death rate of daphnias, depending on their age: delta(age)

#TODO: Plot the best fits for delta for all populations, pchisq() for testing the goodnes of our fit? - ask

#Libraries
library(foreign)   #to read the spss file

#Daphnia cannot live longer than 200 days (in our data not even longer than 123)
lifetime_threshold = 150

#Initialize the variables
age_axis <- c(0:lifetime_threshold)

#Our polynomial Ansatz for delta
delta <- function(parameters, age){
  local_delta <- 0
  if(length(parameters) == 2){
    local_delta <- parameters[1]*age + parameters[2]
  }else if(length(parameters) == 3){
    local_delta <- parameters[1] * age^2 + parameters[2]*age + parameters[3]
  }else if(length(parameters) == 4){
    local_delta <- parameters[1] * age^3 + parameters[2]*age^2 + parameters[3]*age + parameters[4]
  }else if(length(parameters) == 5){
    local_delta <- parameters[1] * age^4 + parameters[2]*age^3 + parameters[3]*age^2 + parameters[4]*age + parameters[5]
  }else if(length(parameters) == 6){
    local_delta <- parameters[1] * age^5 + parameters[2]*age^4 + parameters[3]*age^3 + parameters[4]*age^2 + parameters[5]*age + parameters[6]
  }
 return (local_delta)
}

#Reading the data
dataset <- read.spss("/home/petr/Documents/Master_thesis/AgeEffectsVirulence.sav", to.data.frame = TRUE, use.value.labels = FALSE)
parasyte <- as.character(dataset$Pasteuria)
infecteds <- as.numeric(dataset$Infected)

#Considering control population, infecteds and exposed (but uninfected) separately

control_population <- dataset[parasyte == "Control ",] # Careful, the value of the variable contains the space after Control, possibly FIX the dataset later
host_longevity_control <- control_population$HostLongevity

infected_population <- dataset[infecteds == 1,]
host_longevity_infecteds <-infected_population$HostLongevity

exposed_population <- dataset[(infecteds == 0) & (parasyte != "Control "),]
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

#Likelihood function (log_likelihood) = log(Pruduct over individuals (Product over timesteps(probability of dying in that timestep))) -> log turns pruducts to sums
likelihood_general <- function(parameters, longevity){
  likelihood_local <- 0
  if(min(delta(parameters, age_axis))>=0){
      #does the binomial factor play any role in the likelihood? It does not, right? Since under the logarithm it only plays role of additive constant (independent of delta)
      for (i in 1:length(longevity)) {
        likelihood_local <- likelihood_local -sum(delta(parameters, age_axis[1:longevity[i]])) + log(delta(parameters,longevity[i])) #instead of delta put in the right probability of dying at age longevity[i]
      }
  
    }else{
      likelihood_local <- -Inf
    }
  
  return(-(likelihood_local + log(factorial(length(longevity))))) #optim is minimizing, so we return sign flipped value to get maximization
}


######################## Optimization
## current version of optimization takes ~20 mins
const_initial = c(0.1,0.001,0.00001)
lin_initial = c(0.1,0.001,0.00001)
quadrat_initial = c(0.01,0.0001,0.000001)
cubic_initial = c(0.0001,0.000001,0.00000001)
quartic_initial = c(0.00001,0.0000001,0.000000001)

optimized_likelihood_control_linear <- Inf
optimized_likelihood_infecteds_linear <- Inf
optimized_likelihood_exposed_linear <- Inf

optimized_likelihood_control_quadratic <- Inf
optimized_likelihood_infecteds_quadratic <- Inf
optimized_likelihood_exposed_quadratic <- Inf

optimized_likelihood_control_cubic <- Inf
optimized_likelihood_infecteds_cubic <- Inf
optimized_likelihood_exposed_cubic <- Inf

optimized_likelihood_control_quartic <- Inf
optimized_likelihood_infecteds_quartic <- Inf
optimized_likelihood_exposed_quartic <- Inf

counter <- 0

for (i in const_initial) {
  
  for (j in lin_initial){
    counter <- counter + 1
    message("-------------------------------------------------------------",counter, "/9","-------------------------------------------------------------")
    
    if(optim(par = c(j, i), fn=likelihood_general, longevity = host_longevity_control)$value < optimized_likelihood_control_linear){
      optimized_parameters_control_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = host_longevity_control)$par
      optimized_likelihood_control_linear <- optim(c(j, i), fn=likelihood_general, longevity = host_longevity_control)$value
      message("correction achieved for control population and i=",i," and j=", j)
    }
    
    if(optim(par = c(j, i), fn=likelihood_general, longevity = host_longevity_infecteds)$value < optimized_likelihood_infecteds_linear){
      optimized_parameters_infecteds_linear <- optim(c(j, i), fn=likelihood_general, longevity = host_longevity_infecteds)$par
      optimized_likelihood_infecteds_linear <- optim(c(j, i), fn=likelihood_general, longevity = host_longevity_infecteds)$value
      message("correction achieved for infected population and i=",i," and j=", j)
    }
    
    if(optim(par = c(j, i), fn=likelihood_general, longevity = host_longevity_exposed)$value < optimized_likelihood_exposed_linear){
      optimized_parameters_exposed_linear <- optim(c(j,i), fn=likelihood_general, longevity = host_longevity_exposed)$par
      optimized_likelihood_exposed_linear <- optim(c(j,i), fn=likelihood_general, longevity = host_longevity_exposed)$value
      message("correction achieved for exposed population and i=",i," and j=", j)
    }
    
    for (k in quadrat_initial) {
      if(optim(c(k,j,i), fn=likelihood_general, longevity = host_longevity_control)$value < optimized_likelihood_control_quadratic){
        optimized_parameters_control_quadratic <- optim(c(k,j,i), fn=likelihood_general, longevity = host_longevity_control)$par
        optimized_likelihood_control_quadratic <- optim(c(k,j,i), fn=likelihood_general, longevity = host_longevity_control)$value
        message("correction achieved for control population and i=",i," and j=", j," and k=",k)
      }
      
      if(optim(c(k,j,i), fn=likelihood_general, longevity = host_longevity_infecteds)$value < optimized_likelihood_infecteds_quadratic){
        optimized_parameters_infecteds_quadratic <- optim(c(k,j,i), fn=likelihood_general, longevity = host_longevity_infecteds)$par
        optimized_likelihood_infecteds_quadratic <- optim(c(k,j,i), fn=likelihood_general, longevity = host_longevity_infecteds)$value
        message("correction achieved for infected population and i=",i," and j=", j," and k=",k)
      }
      
      if(optim(c(k,j,i), fn=likelihood_general, longevity = host_longevity_exposed)$value < optimized_likelihood_exposed_quadratic){
        optimized_parameters_exposed_quadratic <- optim(c(k,j,i), fn=likelihood_general, longevity = host_longevity_exposed)$par
        optimized_likelihood_exposed_quadratic <- optim(c(k,j,i), fn=likelihood_general, longevity = host_longevity_exposed)$value
        message("correction achieved for exposed population and i=",i," and j=", j," and k=",k)
      }
      
      for (l in cubic_initial) {
        if(optim(c(l, k, j, i), fn=likelihood_general, longevity = host_longevity_control)$value < optimized_likelihood_control_cubic){
          optimized_parameters_control_cubic <- optim(c(l, k, j, i), fn=likelihood_general, longevity = host_longevity_control)$par
          optimized_likelihood_control_cubic <- optim(c(l, k, j, i), fn=likelihood_general, longevity = host_longevity_control)$value
          message("correction achieved for control population and i=",i," and j=", j," and k=",k," and l=",l)
        }
        
        if(optim(c(l, k, j, i), fn=likelihood_general, longevity = host_longevity_infecteds)$value < optimized_likelihood_infecteds_cubic){
          optimized_parameters_infecteds_cubic <- optim(c(l, k, j, i), fn=likelihood_general, longevity = host_longevity_infecteds)$par
          optimized_likelihood_infecteds_cubic <- optim(c(l, k, j, i), fn=likelihood_general, longevity = host_longevity_infecteds)$value
          message("correction achieved for infected population and i=",i," and j=", j," and k=",k," and l=",l)
        }
        
        if(optim(c(l, k, j, i), fn=likelihood_general, longevity = host_longevity_exposed)$value < optimized_likelihood_exposed_cubic){
          optimized_parameters_exposed_cubic <- optim(c(l, k, j, i), fn=likelihood_general, longevity = host_longevity_exposed)$par
          optimized_likelihood_exposed_cubic <- optim(c(l, k, j, i), fn=likelihood_general, longevity = host_longevity_exposed)$value
          message("correction achieved for exposed population and i=",i," and j=", j," and k=",k," and l=",l)
        }
        for (m in quartic_initial) {
          if(optim(c(m, l, k, j, i), fn=likelihood_general, longevity = host_longevity_control)$value < optimized_likelihood_control_quartic){
            optimized_parameters_control_quartic <- optim(c(m, l, k, j, i), fn=likelihood_general, longevity = host_longevity_control)$par
            optimized_likelihood_control_quartic <-optim(c(m, l, k, j, i), fn=likelihood_general, longevity = host_longevity_control)$value
            message("correction achieved for control population and i=",i," and j=", j," and k=",k," and l=",l," and m=",m)
          }
          
          if(optim(c(m, l, k, j, i), fn=likelihood_general, longevity = host_longevity_infecteds)$value < optimized_likelihood_infecteds_quartic){
            optimized_parameters_infecteds_quartic <- optim(c(m, l, k, j, i), fn=likelihood_general, longevity = host_longevity_infecteds)$par
            optimized_likelihood_infecteds_quartic <-optim(c(m, l, k, j, i), fn=likelihood_general, longevity = host_longevity_infecteds)$value
            message("correction achieved for infected population and i=",i," and j=", j," and k=",k," and l=",l," and m=",m)
          }
          if(optim(c(m, l, k, j, i), fn=likelihood_general, longevity = host_longevity_exposed)$value < optimized_likelihood_exposed_quartic){
            optimized_parameters_exposed_quartic <- optim(c(m, l, k, j, i), fn=likelihood_general, longevity = host_longevity_exposed)$par
            optimized_likelihood_exposed_quartic <- optim(c(m, l, k, j, i), fn=likelihood_general, longevity = host_longevity_exposed)$value
            message("correction achieved for exposed population and i=",i," and j=", j," and k=",k," and l=",l,"a nd m=",m)
          }
        }
      }
    }
  }
}


delta_values_control_linear <- delta(optimized_parameters_control_linear, age_axis)
delta_values_infecteds_linear <- delta(optimized_parameters_infecteds_linear, age_axis)
delta_values_exposed_linear <- delta(optimized_parameters_exposed_linear, age_axis)

delta_values_control_quadratic <- delta(optimized_parameters_control_quadratic, age_axis)
delta_values_infecteds_quadratic <- delta(optimized_parameters_infecteds_quadratic, age_axis)
delta_values_exposed_quadratic <- delta(optimized_parameters_exposed_quadratic, age_axis)

delta_values_control_cubic <- delta(optimized_parameters_control_cubic, age_axis)
delta_values_infecteds_cubic <- delta(optimized_parameters_infecteds_cubic, age_axis)
delta_values_exposed_cubic <- delta(optimized_parameters_exposed_cubic, age_axis)

delta_values_control_quartic <- delta(optimized_parameters_control_quartic, age_axis)
delta_values_infecteds_quartic <- delta(optimized_parameters_infecteds_quartic, age_axis)
delta_values_exposed_quartic <- delta(optimized_parameters_exposed_quartic, age_axis)

############ TRY 5th order polynomial
# fifth_initial <- c(0.0000001,0.000000001,0.0000000001)
# optimized_likelihood_control_fifth_order <- Inf
# for (i in const_initial) {
#   
#   for (j in lin_initial){
#     
#     for (k in quadrat_initial) {
#       
#       
#       for (l in cubic_initial) {
#         
#         for (m in quartic_initial) {
#           for (fifth in fifth_initial) {
#             counter <- counter + 1
#             message("-------------------------------------------------------------",counter, "/729","-------------------------------------------------------------")
#             if(optim(par=c(fifth,optimized_parameters_control_quartic), fn=likelihood_general, longevity = host_longevity_control)$value < optimized_likelihood_control_fifth){
#               optimized_parameters_control_fifth_order <- optim(par=c(fifth,optimized_parameters_control_quartic), fn=likelihood_general, longevity = host_longevity_control)$par
#               optimized_likelihood_control_fifth_order <- optim(par=c(fifth,optimized_parameters_control_quartic), fn=likelihood_general, longevity = host_longevity_control)$value
#             }
#           }
#         }
#       }
#     }
#   }
# }


#check overlap with survival data

survival_model <- function(delta, x = age_axis){
  exponent <- c(rep(0,length(x)))
  for (a in x) {
    exponent[a+1] <- -sum(delta[1:(a+1)])
  }
  return(exp(exponent))
}

#########################   Outputs and plotting

#parameters for delta, after optimization

optimized_parameters_control_quartic
optimized_parameters_infecteds_quartic
optimized_parameters_exposed_quartic

optimized_parameters_control_cubic
optimized_parameters_infecteds_cubic
optimized_parameters_exposed_cubic

optimized_parameters_control_quadratic
optimized_parameters_infecteds_quadratic
optimized_parameters_exposed_quadratic

optimized_parameters_control_linear
optimized_parameters_infecteds_linear
optimized_parameters_exposed_linear

########Likelihood Ratio Test
LRT_control_lin_to_quadratic <- 2*(optimized_likelihood_control_linear - optimized_likelihood_control_quadratic)
LRT_control_quadratic_to_cubic <- 2*(optimized_likelihood_control_quadratic - optimized_likelihood_control_cubic)
LRT_control_cubic_to_quartic <- 2*(optimized_likelihood_control_cubic - optimized_likelihood_control_quartic)
#LRT_control_quartic_to_fifth_order <- 2*(optimized_likelihood_control_quartic - optimized_likelihood_control_fifth_order)

LRT_control_lin_to_quadratic
LRT_control_quadratic_to_cubic
LRT_control_cubic_to_quartic
#LRT_control_quartic_to_fifth_order

LRT_infecteds_lin_to_quadratic <- 2*(optimized_likelihood_infecteds_linear - optimized_likelihood_infecteds_quadratic)
LRT_infecteds_quadratic_to_cubic <- 2*(optimized_likelihood_infecteds_quadratic - optimized_likelihood_infecteds_cubic)
LRT_infecteds_cubic_to_quartic <- 2*(optimized_likelihood_infecteds_cubic - optimized_likelihood_infecteds_quartic)
LRT_infecteds_quadratic_to_quartic <- 2*(optimized_likelihood_infecteds_quadratic -optimized_likelihood_infecteds_quartic)


LRT_infecteds_lin_to_quadratic
LRT_infecteds_quadratic_to_cubic
LRT_infecteds_cubic_to_quartic
LRT_infecteds_quadratic_to_quartic

LRT_exposed_lin_to_quadratic <- 2*(optimized_likelihood_exposed_linear - optimized_likelihood_exposed_quadratic)
LRT_exposed_quadratic_to_cubic <- 2*(optimized_likelihood_exposed_quadratic - optimized_likelihood_exposed_cubic)
LRT_exposed_cubic_to_quartic <- 2*(optimized_likelihood_exposed_cubic - optimized_likelihood_exposed_quartic)

LRT_exposed_lin_to_quadratic
LRT_exposed_quadratic_to_cubic
LRT_exposed_cubic_to_quartic

#pchisq() test of the LR for control population
pchisq(LRT_control_lin_to_quadratic,df = 1,lower.tail = FALSE)
pchisq(LRT_control_quadratic_to_cubic,df = 1,lower.tail = FALSE)
pchisq(LRT_control_cubic_to_quartic,df = 1,lower.tail = FALSE)
#pchisq(LRT_control_cubic_to_fifth_order, df = 1, lower.tail = FALSE) #likelihood is insufficient


#pchisq() test of the LR for infected population
pchisq(LRT_infecteds_lin_to_quadratic,df = 1,lower.tail = FALSE)
pchisq(LRT_infecteds_quadratic_to_cubic,df = 1,lower.tail = FALSE)
pchisq(LRT_infecteds_cubic_to_quartic,df = 1,lower.tail = FALSE)
pchisq(LRT_infecteds_quadratic_to_quartic, df = 2, lower.tail = FALSE)

#pchisq() test of the LR for exposed population
pchisq(LRT_exposed_lin_to_quadratic,df = 1,lower.tail = FALSE)
pchisq(LRT_exposed_quadratic_to_cubic,df = 1,lower.tail = FALSE)
pchisq(LRT_exposed_cubic_to_quartic,df = 1,lower.tail = FALSE)


#plotting the survival compared to linear and quadratic model of delta
pdf("Control_population_survival.pdf")
plot_control <- plot(age_axis, survivals_control, col="blue", xlab = "Age [days]", ylab = "Survivals in control population [fraction]")
lines(age_axis, survival_model(delta = delta_values_control_linear), col="red")
lines(age_axis, survival_model(delta = delta_values_control_quadratic), col="green")
lines(age_axis, survival_model(delta = delta_values_control_cubic), col="yellow")
lines(age_axis, survival_model(delta = delta_values_control_quartic), col="purple")
#lines(age_axis, survival_model(delta = delta(optimized_parameters_control_fifth_order, age = age_axis)), col="red")
legend(0,0.3, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("blue", "red", "green", "yellow", "purple"), lty = 1, cex = 0.8)
dev.off()

pdf("Infected_population_survival.pdf")
plot_infecteds <- plot(age_axis, survivals_infecteds, col="blue", xlab = "Age [days]", ylab = "Survivals in infected population [fraction]")
lines(age_axis, survival_model(delta = delta_values_infecteds_linear), col="red")
lines(age_axis, survival_model(delta = delta_values_infecteds_quadratic), col="green")
lines(age_axis, survival_model(delta = delta_values_infecteds_cubic), col="yellow")
lines(age_axis, survival_model(delta = delta_values_infecteds_quartic), col="purple")
legend(60,1, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("blue", "red", "green", "yellow", "purple"), lty = 1, cex = 0.8)
dev.off()


pdf("Exposed_population_survival.pdf")
plot_exposed <- plot(age_axis, survivals_exposed, col="blue", xlab = "Age [days]", ylab = "Survivals in exposed population [fraction]")
lines(age_axis, survival_model(delta = delta_values_exposed_linear), col="red")
lines(age_axis, survival_model(delta = delta_values_exposed_quadratic), col="green")
lines(age_axis, survival_model(delta = delta_values_exposed_cubic), col="yellow")
lines(age_axis, survival_model(delta = delta_values_exposed_quartic), col="purple")
legend(0,0.3, legend = c("Data", 'Linear model', "Quadratic model", "Cubic model", "4th order polynomial model"), col = c("blue", "red", "green", "yellow", "purple"), lty = 1, cex = 0.8)
dev.off()

optimized_parameters_control_quadratic
delta(optimized_parameters_control_quadratic, age_axis)

#plotting delta, for different populations

pdf("Quadratic_models.pdf")
plot_deltas_quadratic <- plot(age_axis, delta(optimized_parameters_control_quadratic, age_axis), xlab = "Age [days]", ylab = "Delta under quadratic approx.", col = "blue", type = "l")
lines(age_axis, delta(optimized_parameters_infecteds_quadratic, age_axis), col = "red")
lines(age_axis, delta(optimized_parameters_exposed_quadratic, age_axis), col = "green")
legend(100,0.02, legend = c("Control", 'Infected', "Exposed"), col = c("blue", "red", "green"), lty = 1)
dev.off()

pdf("Linear_models.pdf")
plot_deltas_linear <- plot(age_axis, delta(optimized_parameters_control_linear, age_axis), xlab = "Age [days]", ylab = "Delta under linear approx.", col = "blue", type = "l")
lines(age_axis, delta(optimized_parameters_infecteds_linear, age_axis), col = "red")
lines(age_axis, delta(optimized_parameters_exposed_linear, age_axis), col = "green")
legend(100,0.01, legend = c("Control", 'Infected', "Exposed"), col = c("blue", "red", "green"), lty = 1)
dev.off()

pdf("Cubic_models.pdf")
plot_deltas_cubic <- plot(age_axis, delta(optimized_parameters_control_cubic, age_axis), xlab = "Age [days]", ylab = "Delta under cubic approx.", col = "blue", type = "l")
lines(age_axis, delta(optimized_parameters_infecteds_cubic, age_axis), col = "red")
lines(age_axis, delta(optimized_parameters_exposed_cubic, age_axis), col = "green")
legend(0,0.30, legend = c("Control", 'Infected', "Exposed"), col = c("blue", "red", "green"), lty = 1)
dev.off()

pdf("Quartic_models.pdf")
plot_deltas_quartic <- plot(age_axis, delta(optimized_parameters_control_quartic, age_axis), xlab = "Age [days]", ylab = "Delta under quartic approx.", col = "blue", type = "l")
lines(age_axis, delta(optimized_parameters_infecteds_quartic, age_axis), col = "red")
lines(age_axis, delta(optimized_parameters_exposed_quartic, age_axis), col = "green")
legend(0,0.30, legend = c("Control", 'Infected', "Exposed"), col = c("blue", "red", "green"), lty = 1)
dev.off()

pdf("Chosen_models.pdf")
plot_deltas <- plot(age_axis, delta(optimized_parameters_control_quartic, age_axis), xlab = "Age [days]", ylab = "Delta for various populations", col = "blue", type = "l")
lines(age_axis, delta(optimized_parameters_infecteds_quadratic, age_axis), col = "red")
lines(age_axis, delta(optimized_parameters_exposed_quartic, age_axis), col = "green")
lines(age_axis, delta(optimized_parameters_infecteds_quartic, age_axis), col = "purple")
legend(0,0.40, legend = c("Control (quartic)", 'Infected(quadratic)', "Exposed(quartic)", "Infected(quartic)"), col = c("blue", "red", "green", "purple"), lty = 1, cex = 0.8)
dev.off()

pdf("Virulence_97days.pdf")
plot_virulence <- plot(age_axis[1:97], delta(optimized_parameters_infecteds_quartic,age_axis[1:97])-delta(optimized_parameters_control_quartic, age_axis[1:97]), xlab = "Age [days]", ylab = "Virulence for 4th order fit of deathrate of infecteds", col = "blue", type = "l")
lines(age_axis[1:97], delta(optimized_parameters_infecteds_quadratic, age_axis[1:97])-delta(optimized_parameters_control_quartic, age_axis[1:97]), col = "red")
legend(5,0.13, legend = c("Virulence (quartic/quartic)", 'Virulence(quartic/quadratic)'), col = c("blue", "red"), lty = 1, cex = 0.8)
dev.off()
View(dataset)