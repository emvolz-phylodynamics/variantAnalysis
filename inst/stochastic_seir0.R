# stochastic SEIR 
update(cI)  <- cI + n_SI
update(E) <- E + n_SI  - n_EI
update(I) <- I + n_EI - n_IR 

## flux
dx <- 1

En_SI <- R*gamma2*I*dx
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
R <- user(1.2) 
gamma1 <- user( 1/4 ) ## TODO update defaults
gamma2 <- user( 1/3 )
