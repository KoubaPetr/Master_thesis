#new script for testing of the implementation of the data into tge age-model
max_age <- 199
age.model <-
  function(D=c(1,rep(0,max_age)),
           I=c( 1,rep(0,max_age)),
           P=0, #should we simulate the transmission after death by pathogen presence in the environment?
           t.end=150){
    
    ## carrying capacity:
    K <- 1000
    
    ## D is the daphnia vector, one entry per age group;
    #D <- vector("numeric", length=200)
    ## we define age groups for 200 days because daphnia can live half a year in captivity
    
    ## I is the infected daphnia vector, one entry per age group;
    #I <- vector("numeric", length=200)
    
    # delta function corresponding to the fit for the control population in the experiment
    
    delta_testing <- function(a){
      return(7.360907e-09*a^4 -1.134501e-06*a^3 + 4.977512e-05*a^2 -4.648003e-04*a + 1.241934e-03)
    }
       
    delta_testing_infected_quadratic <- function(a){
      return(2.467349e-05*a^2 -6.657129e-04*a + 4.484470e-03)
    }
    
    delta_testing_infected_quartic<-function(a){
      return(-4.608850e-09*a^4 + 8.041395e-07*a^3 -1.759914e-05*a^2 + 1.014732e-04*a + 2.261560e-05)
    }
    
    plot(0:100, delta_testing_infected_quartic(0:100)-delta_testing(0:100)) #plotting virulence for the first 101 days
    ## age dependent reproduction:
    r <- c(#rep(0,10),
           #3,4,5,5,5,5,5,5,5,5,
           rep(0,200))
    ## age dependent mortaility: computed based on the fit for control population
    m <- delta_testing(a = 0:max_age)
    m
    ## age dependent susceptibility:
    f <- c(#rep(0.1, 5),
           #rep(0.05, 10),
           #rep(0.02, 15),
           #rep(0.02,170),
          rep(0,200))
    ## infectiousness of a burst (ie hte number of parasite spores that come out
    ## of an infected flea at death):
    b <- 0.1
    ## age dependent virulence: (there was no difference between the ages in virulence, was there?)
    v <- delta_testing_infected_quadratic(a = 0:max_age) -m
    
    #    D[1] <- 10 # starting cohort
    #    I[10] <- 2 # starting cohort infecteds
    #    P <- 0
    output <- data.frame(t(c(0,P,D,I)))
    names(output) <- c("day","P",
                       paste0("D",1:length(D)),
                       paste0("I",1:length(I)))
    
    for(t in 1:t.end){
      ##cat(paste0("t=",t,"\n")) # for diagnostics
      ## let fleas reproduce:
      newborns <- sum(r*D)
      ## dying:
      D <- (1-m)*D
      ## age every flea:
      D.new <- c(newborns,D[-length(D)]) ####missing loss of susceptibles due to infection?
      
      ## dying by infection: ##### and also mortality on the background
      I <- (1-v-m)*I
      P <- sum(v*I) ##### rewriting the pathogens, rather incrementing, right? Or in this simple model do the pathogens live just one day in the environment?
      ## new infections:
      I.new <- b*P*f*D
      ## age every flea:
      I <- I.new + c(0,I[-length(I)])
      
      ## shrink to carrying capacity:
      D <- D.new
      if(sum(I)+sum(D) > K){
        D <- D * K/(sum(I)+sum(D))
        I <- I * K/(sum(I)+sum(D))
      }
      output <- rbind(output,
                      c(t,P,D,I))
    }
    output
  }
am <- age.model()

survival_model <- function(delta, x){
  exponent <- c(rep(0,length(x)))
  for (a in x) {
    exponent[a+1] <- -sum(delta[1:(a+1)])
  }
  return(exp(exponent))
}

plotting_function_per_day <- 
    function(model,day,population){
      indexing <- 0
      if(population == "D"){
        indexing <- 2
      }else if(population == "I"){
        indexing <- 2 + max_age
      }
      plot(0:150,t(unname(as.vector(model[indexing + day]))))
    }
plotting_function_per_group <-
  function(model = am, population = "D"){
    indexing <- 0
    survivals_local <- c(rep(0,max_age+1))
    
    if(population == "D"){
      indexing <- 2
    }else if(population == "I"){
      indexing <- 2 + max_age + 1
    }
    
    for (i in 1:length(model$D1)) {
      survivals_local[i] <- t(model[i+indexing])[i]
    }
    #survivals_local
    if(population == "D"){
      age_axis_MLE <- 0:max_age
      delta_values_control_MLE <-delta_testing(age_axis_MLE)
      pdf("Comparison_of_models_control.pdf")
      plot(0:max_age, survivals_local,xlab = "Age a = Time t [days]", ylab = "Survivals in control population [fraction]")
      lines(age_axis_MLE, survival_model(x = age_axis_MLE, delta = delta_values_control_MLE), col="purple", lwd = 3)
      legend(47,1, legend = c("Model from Section 2", 'Survival using the optimal delta'), col = c("black", "purple"), lwd =3, lty = 1, cex = 0.8)
      dev.off()
    }else if(population == "I"){
      age_axis_MLE <- 0:max_age
      delta_values_control_MLE <-delta_testing_infected_quadratic(age_axis_MLE)
      pdf("Comparison_of_models_infecteds.pdf")
      plot(0:max_age, survivals_local,xlab = "Age a = Time t [days]", ylab = "Survivals in infected population [fraction]")
      lines(age_axis_MLE, survival_model(x = age_axis_MLE, delta = delta_values_control_MLE), col="purple", lwd = 3)
      legend(47,1, legend = c("Model from Section 2", 'Survival using the optimal delta'), col = c("black", "purple"), lwd =3, lty = 1, cex = 0.8)
      dev.off()
    }
  }
# am$D145[145] #why this value??
# t(am[147])[145]

#lines(age_axis, survival_model(delta = delta_values_control_quartic), col="purple")

plotting_function_per_day(model = am, day = 20,population = "D")
plotting_function_per_group()
plotting_function_per_group(model = am, population = "I")