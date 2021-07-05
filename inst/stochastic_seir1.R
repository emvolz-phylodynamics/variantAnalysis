# stochastic SEIR 
update(cI)  <- cI + n_SI
update(E) <- E + n_SI  - n_EI
update(I) <- I + n_EI - n_IR 

## flux
dx <- 1
pini <- max(0, min(1, (step-tini)/(tequil-tini) )) #TODO does not work
Rt <- Requil * pini  + R * (1-pini)

En_SI <- Rt*gamma2*I*dx
n_SI <- max(0, rnorm(En_SI, sqrt(En_SI) )  )

En_EI <- gamma1 * E * dx 
n_EI <- max(0, rnorm( En_EI, sqrt(En_EI )))

En_IR <- gamma2 * I * dx 
n_IR <- max(0, rnorm( En_IR, sqrt( En_IR )))


## initial 
initial(cI) <- 1
initial(E) <- 1
initial(I) <- 0 

## inputs
R <- user(1.75) 
Requil <- user( 1.25 )
gamma1 <- user( 1/5 ) ## matches SEIJR models
gamma2 <- user( 1/3 )
tini <- user(0)
tequil <- user( 45 )
