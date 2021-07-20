# v2 incl super-spreading 

# stochastic SEIR 
update(cI)  <- cI + n_SE
update(E) <- E + n_SE  - n_EIl - n_EIh 
update(Il) <- Il + n_EIl - n_IlR 
update(Ih) <- Ih + n_EIh - n_IhR 

## flux
dx <- 1

## R changes over time 
pini <- max(0, min(1, (step-tini)/(tequil-tini) )) #TODO does not work
Rt <- Requil * pini  + R * (1-pini)

## different R in high and low risk pops
Rl <- Rt / ( ph*tau + 1 - ph )
Rh <- Rl * tau 

## transmission 
En_SE <- (Rh*gamma2*Ih + Rl*gamma2*Il)*dx
n_SE <- max(0, rnorm(En_SE, sqrt(En_SE) )  )

En_EI <- gamma1 * E * dx 
n_EI <- max(0, rnorm( En_EI, sqrt(En_EI )))

## break down to Il and Ih 
### random proportion going to h 
phx <- max(0, min(1 , rnorm( ph , sqrt( ph*(1-ph)/n_EI )  ) ))
n_EIl <- (1-phx) * n_EI
n_EIh <- phx * n_EI

## I -> R
En_IlR <- gamma2 * Il * dx 
n_IlR <- max(0, rnorm( En_IlR, sqrt( En_IlR )))

En_IhR <- gamma2 * Ih * dx 
n_IhR <- max(0, rnorm( En_IhR, sqrt( En_IhR )))


## initial , approx distribution if Rt is close to one
initial(cI) <- EI0
#~ initial(E) <- 1
#~ initial(Ih) <- 0 
#~ initial(Il) <- 0 
initial(E) <- EI0 * (1/gamma1) / (1/gamma1 + 1/gamma2) 
initial(Ih) <- EI0 * ph * (1/gamma2) / (1/gamma1 + 1/gamma2 )
initial(Il) <- EI0 * (1-ph) * (1/gamma2) / (1/gamma1 + 1/gamma2 )

## inputs
R <- user( 1.75 ) 
Requil <- user( 1.25 ) 
gamma1 <- user( 1/5 )  ## matches SEIJR models
gamma2 <- user( 1/3 ) 
tini <- user( 0 ) 
tequil <- user( 45 ) 
ph <- user( .2 ) 
tau <- user( 74 ) 
EI0 <- user(1)   # initial infections
