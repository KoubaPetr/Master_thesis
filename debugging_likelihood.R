#Debugging likelihoog

require(foreign)   #to read the spss file

#Daphnia cannot live longer than 200 days (in our data not even longer than 123)
lifetime_threshold = 150

#Initialize the variables ---------------------------------------------------------------
age_axis <- c(0:lifetime_threshold)

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

delta_weibull <- function(parameters, age){
  k <- parameters[1]
  lambda <- parameters[2]
  return(k/lambda * (age/lambda)^(k-1))
  
}

#Reading the data -----------------------------------------------------------------------
dataset <- read.spss("/home/petr/Documents/Master_thesis/AgeEffectsVirulence.sav", to.data.frame = TRUE, use.value.labels = FALSE)
parasyte <- as.character(dataset$Pasteuria)
infecteds <- as.numeric(dataset$Infected)
aai <- as.numeric(dataset$InfectionAge)

exposed_population <- dataset[(infecteds == 0) & (parasyte != "Control "),]
longevity_exposed <- as.numeric(exposed_population$HostLongevity)
survivals_exposed <- rep(0,lifetime_threshold+1)

for (x in age_axis) {
  survivals_exposed[x+1] <- sum(longevity_exposed>x)/length(longevity_exposed)
}

### log likelihood TODO: Rewrite the second cycle analytically
# likelihood_subroutine <- function(par, d, a){
#   return(- sum(1:(d-a) * delta(par, (a-1):(d-1))) + log(delta(par,d)))
# }

 likelihood_general <- function(p, longevity, age_at_infection){
   likelihood_local <- 0
   if(min(delta(parameters = p, age_axis))>=0){
     for (d in longevity) {
       likelihood_local <- likelihood_local -sum(delta(parameters = p,(age_at_infection-1):d)) + log(delta(parameters = p, d))
       #likelihood_local <- likelihood_local - sum(1:(d-age_at_infection) * delta(parameters, (age_at_infection-1):(d-1))) + log(delta(parameters,d))
       #sum(sapply(longevity, likelihood_subroutine, par = parameters, a = age_at_infection))#}
     }
   }else{
     likelihood_local <- -Inf
   }

   return(-likelihood_local) #optim is minimizing, so we return sign flipped value to get maximization #consider adding the binomial factor to have absolute likelihood,
   #keep in mind the different ages of infections which play role
 }



start <- Sys.time()
#### Optimization
optimized_likelihood_exposed_quartic <- Inf

const_initial = c(0.1,0.001,0.00001)
#const_initial = c(0.00001)
lin_initial = c(0.1,0.001,0.00001)
#lin_initial = c(0.00001)
quadrat_initial = c(0.01,0.0001,0.000001)
#quadrat_initial = c(0.000001)
cubic_initial = c(0.0001,0.000001,0.00000001)
#cubic_initial = c(0.00000001)
quartic_initial = c(0.00001,0.0000001,0.000000001)
#quartic_initial = c(0.000000001)


# # weibull optimization
# optimized_likelihood_exposed_weib <- Inf
# 
# for (i in 1:9) {
#   for (j in 1:9) {
#     
#     if(optim(par = c(10^(i-1),10^(j-1)), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$value < optimized_likelihood_exposed_weib){
#             optimized_likelihood_exposed_weib <- optim(par = c(10^(i-1),10^(j-1)), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$value
#             optimized_parameters_exposed_weib <- optim(par = c(10^(i-1),10^(j-1)), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$par
#             message("correction achieved for uninfected population and k= ",i," and lambda= ", j)
#           }
#   }
# }



 for (i in const_initial) {

   for (j in lin_initial){

     # if(optim(par = c(j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$value < optimized_likelihood_exposed_linear){
     #   optimized_likelihood_exposed_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$value
     #   optimized_parameters_exposed_linear <- optim(par = c(j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$par
     #   message("correction achieved for uninfected population and i=",i," and j=", j)
     # }

     for (k in quadrat_initial) {

       # if(optim(par = c(k,j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$value < optimized_likelihood_exposed_quadratic){
       #   optimized_likelihood_exposed_quadratic <- optim(par = c(k,j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$value
       #   optimized_parameters_exposed_quadratic <- optim(par = c(k,j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$par
       #   message("correction achieved for uninfected population and i=",i," and j=", j, "and k=",k)
       # }
       for (l in cubic_initial) {

         # if(optim(par = c(l,k,j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$value < optimized_likelihood_exposed_cubic){
         #   optimized_likelihood_exposed_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$value
         #   optimized_parameters_exposed_cubic <- optim(par = c(l, k, j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$par
         #   message("correction achieved for uninfected population and i=",i," and j=", j, "and k=",k," and l=",l)
         # }

         for (m in quartic_initial) {

           if(optim(par = c(m,l,k,j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$value < optimized_likelihood_exposed_quartic){
             optimized_likelihood_exposed_quartic <- optim(par = c(m,l, k, j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$value
             optimized_parameters_exposed_quartic <- optim(par = c(m,l, k, j, i), fn=likelihood_general, longevity = longevity_exposed, age_at_infection = 1)$par
             message("correction achieved for uninfected population and i=",i," and j=", j, "and k=",k," and l=",l, " and m=",m)
           }
         }
       }
     }
   }
 }

#### Survival model
survival_model <- function(delta, x = age_axis){
  exponent <- c(rep(0,length(x)))
  for (a in x) {
    exponent[a+1] <- -sum(delta[1:(a+1)])
  }
  return(exp(exponent))
}
stop <- Sys.time()
stop - start
optimized_parameters_exposed_quartic
#### PLotting

#likelihood_fitting_since_age_16 <- optimized_likelihood_exposed_quartic
likelihood_fitting_since_age_16
likelihood_fitting_since_age_1

plot(x = age_axis, y = survivals_exposed, xlab = "Age [Days]", ylab = "Survival_Exposed", type = "s", lwd =2)
lines(x = age_axis, y = survival_model(delta = delta(parameters = optimized_parameters_exposed_quartic, age = age_axis), x = age_axis), col = "red")
legend(90,1, legend = c("Data", "4th order polynomial model"), col = c("black", "red"), lty = 1, cex = 0.8)
text(50, 0.4, "version fitting since age 1, alive till d")

# plot(x = age_axis, y = survivals_exposed, xlab = "Age [Days]", ylab = "Survival_Exposed", type = "s", lwd =2)
# lines(x = age_axis, y = survival_model(delta = delta_weibull(parameters = optimized_parameters_exposed_weib, age = age_axis), x = age_axis), col = "red")
# legend(90,1, legend = c("Data", "weibull model"), col = c("black", "red"), lty = 1, cex = 0.8)
# plot(delta_weibull(optimized_parameters_exposed_weib, age = age_axis))

#TODO:
# weibull hazard is polynomial, so it is probably useless, we can probably stick with the current delta
# implement the current likelihood into the main code