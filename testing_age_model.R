#new script for testing of the implementation of the data into tge age-model

age.model <-
  function(D=c(10,rep(0,199)),
           I=c( 2,rep(0,199)),
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
    
    delta <- function(a){
      return(4.420315e-07*a^4 - 7.621028e-05*a^3 + 3.844089e-03*a^2 - 4.674704e-02*a + 1.944348e-01)
    }
    
    ## age dependent reproduction:
    r <- c(rep(0,10),
           3,4,5,5,5,5,5,5,5,5,
           rep(0,180))
    ## age dependent mortaility: computed based on the fit for control population
    m <- delta(a = 0:199)
    
    ## age dependent susceptibility:
    f <- c(rep(0.1, 5),
           rep(0.05, 10),
           rep(0.02, 15),
           rep(0.02,170))
    ## infectiousness of a burst (ie hte number of parasite spores that come out
    ## of an infected flea at death):
    b <- 0.1
    ## age dependent virulence: (there was no difference between the ages in virulence, was there?)
    v <- c(rep(0.05, 200))
    
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
      D.new <- c(newborns,D[-length(D)])
      
      ## dying by infection: #and also mortality on the background
      I <- (1-v-m)*I
      P <- sum(v*I)
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