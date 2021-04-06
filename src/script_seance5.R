## chargement packages
library(deSolve)
library(ggplot2)

##  fonctin SIR de base
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

### parametres
## Initialisation Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S = 1-1e-6, I = 1e-6, R = 0.0)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 0.8387, gamma = 0.1346)
1/0.1346
1/15
## Time frame
times      <- seq(0, 100, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)

plot(times, out$S,type="l",col="blue",ylim=c(0,1),ylab="% SIR")
legend("topright", col=c("blue","red","green"), legend=c("S", "I", "R"), lwd=1)
points(times, out$I,type="l",col="red")
points(times, out$R,type="l",col="green")

# ### diminution du taux de guérison 
# ### parametres
# ## Initialisation Susceptible 0.999999, Infected 0.000001, Recovered 0
# init       <- c(S = 1-1e-6, I = 1e-6, R = 0.0)
# ## beta: infection parameter; gamma: recovery parameter
# parameters <- c(beta = 0.2387, gamma = 1/15)
# # 1/0.1346 7 jours environ 
# 1/15
# ## Time frame
# times      <- seq(0, 400, by = 1)
# 
# ## Solve using ode (General Solver for Ordinary Differential Equations)
# out <- ode(y = init, times = times, func = sir, parms = parameters)
# ## change to data frame
# out <- as.data.frame(out)
# 
# #plot(times, out$S,type="l",col="blue",ylim=c(0,1),ylab="% SIR")
# #legend("topright", col=c("blue","red","green"), legend=c("S", "I", "R"), lwd=1)
# points(times, out$S,type="l",col="blue",lwd=2,lty=2)
# points(times, out$I,type="l",col="red",lwd=2,lty=2)
# points(times, out$R,type="l",col="green",lwd=2,lty=2)
# 
# 
# ### augmentation du taux de contagion 
# ### parametres
# ## Initialisation Susceptible 0.999999, Infected 0.000001, Recovered 0
# init       <- c(S = 1-1e-6, I = 1e-6, R = 0.0)
# ## beta: infection parameter; gamma: recovery parameter
# parameters <- c(beta = 1.2387, gamma = 0.1346)
# 1/0.1346
# ## Time frame
# times      <- seq(0, 400, by = 1)
# 
# ## Solve using ode (General Solver for Ordinary Differential Equations)
# out <- ode(y = init, times = times, func = sir, parms = parameters)
# ## change to data frame
# out <- as.data.frame(out)
# 
# #plot(times, out$S,type="l",col="blue",ylim=c(0,1),ylab="% SIR")
# #legend("topright", col=c("blue","red","green"), legend=c("S", "I", "R"), lwd=1)
# points(times, out$S,type="l",col="blue",lwd=2,lty=2)
# points(times, out$I,type="l",col="red",lwd=2,lty=2)
# points(times, out$R,type="l",col="green",lwd=2,lty=2)


### PRISE EN COMPTE DE LA VACCINATION
##  fonctin SIR avec vaccination
sir_with_vacc <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I              - v * S
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I + v * S
    
    return(list(c(dS, dI, dR)))
  })
}

### parametres
## Initialisation Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S = 1-1e-6, I = 1e-6, R = 0.0)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 0.8387, gamma = 0.1346, v=0.005)
1/0.1346
1/15
## Time frame
times      <- seq(0, 100, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir_with_vacc, parms = parameters)
## change to data frame
out <- as.data.frame(out)

points(times, out$S,type="l",col="blue",lwd=2,lty=2)
points(times, out$I,type="l",col="red",lwd=2,lty=2)
points(times, out$R,type="l",col="green",lwd=2,lty=2)

### récupérer les données covid 19
datacovid <- read.csv2("covid19-du-2020-03-04-au-2020-03-25.csv",sep=",")
max(datacovid[,-1])/70000000
plot(datacovid$Bretagne, ylab="Nombre de cas", xlab="Temps",col=2)
for (i in 3:21) {
  points(datacovid[,i],col=i)
}

## on est dans les 20 premiers jours de l'épidémie 
## S vaut 1 au maximum 
## S vaut
1-max(datacovid[,-1])/70000000
## pour les 20 premiers jours, on peut remplacer S par 1
## le système 
# dS <- -beta * S * I
# dI <-  beta * S * I - gamma * I
#se simplifie pour donner
# dI <-  beta * 1 * I - gamma * I
# la solution
## I(t)=I_o exp((beta-gamma)t)
## pour rappel, au début de l'épidémie Ro= beta/gamma
## I(t)=I_o exp(gamma(beta/gamma-1)t)
## I(t)=I_o exp(gamma(R0-1)t)
## je vois exp donc on a envie de loguer
## log(I(t))= log(I_o)+ gamma(R0-1)t
Idata <- datacovid$Bretagne
vectemps <- 1:20 
## idée : reconnaitre une equation du type
## y=ax+b
## y : log(I(t))
## x : t 
## a : slope= gamma(R0-1)
## b : log(I_o)
reslm <- lm(log(Idata) ~ vectemps)
slope <- coefficients(reslm)[2]
## slope = gamma(R0-1)
## R0 = slope/gamma +1
## gamma =1/temps guérison
gamma=1/15
R0=slope/gamma +1
R0
