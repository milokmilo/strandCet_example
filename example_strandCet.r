# ###################################################################################### #
#  strandCet: R package for estimating natural and anthropogenic mortality-at-age 
#             of cetaceans from age-structured strandings.
#                                   Camilo Saavedra
#           Instituto Español de Oceanografía, Centro Oceanográfico de Vigo
#                           P.O. Box 1552, 36200 Vigo, Spain.
#                       E-mail: camilo.saavedra.penas@gmail.com
# ###################################################################################### #

# Below you can find the commands used to run the example of the paper.
# Note that model fits and estimated parameters could slightly from one run to another 
# General stepwise instructions for the different functions and their options, 
# as well as other examples can be found in the *strandCet* manual
# <https://cran.r-project.org/web/packages/strandCet/strandCet.pdf>, 
# For further package details you can also consult the help files -> help("strandCet").

install.packages("strandCet")
library(strandCet)
install.packages("ggplot2")
library(ggplot2)


#### DATA ####

## Data from Stolen and Barlow (2003) 

# Number of dead dolphins
N <- c(36,14,32,15,13,8,4,5,3,3,2,3,3,6,3,8,5,5,8,3,11,4,1,4,3,3,6,2,1,2,0,1,0,1,0,2)
# Age range of the dead dolphins 
age <- 0:35
    
# Add normally distributed bycaught dolphins with mode before maturity
set.seed(123)
Nbyc <- N + c(sort(round(rnorm(6, mean = 10, sd = 5))),
            rev(sort(round(rnorm(6, mean = 10, sd = 5)))), rep(0,24))

# Create the two datasets
dataSI <- data.frame(age, N) # Original data from Stolen and Barlow (2003) 
dataHP <- data.frame(age, Nbyc) # Theoretical bycaught dolphins added 


#### SILER MODEL ####

## Total mortality-at-age can be calculated using the Siler model
## The Siler model is applied to the original data without bycatch

# Observed life table with the original data
lifeSI <- life.tab(dataSI)

# Siler model with the original data
modSI <- Si.mod(dataSI)

# Predict the Siler model with the original data
predSI <- Si.pred(dataSI, modSI)

# Life table of the Siler model with the original data
life.Siler <- Est.life.tab(Est.qx = predSI$qx.tot, age = 0:35, n = 1000)


#### HELIGMAN-POLLARD MODEL ####

## Total mortality-at-age can be separated in natural and anthropogenic mortality 
## The Heligman-Pollard model is applied to the data with bycatch

# Observed life table with the original data
lifeHP <- life.tab(dataHP)

# Select non-informative priors for the Heligman-Pollard model
priors <- data.frame(priors.lo = c(0,0,0,0,0,0,1,0,1),
                     priors.hi = c(1,5,1,0.01,0.5,10,15,0.01,1.5))

# Compile priors in the required format
q0 <- HP.priors(pri.lo = priors$priors.lo,
                pri.hi = priors$priors.hi,
                theta.dim = 9)

# Run the Heligman-Pollard model
modHP <- HP.mod(prior = q0, lifeTab = lifeHP,
                 K = 10, d = 10, B = 500, CI = 90)

# Predict the Heligman-Pollard model model
predHP <- HP.pred(life = lifeHP, HPout = modHP, age = age)

# Heligman-Pollard plot
HPcurves <- ggplot(predHP, aes(age, qx.tot)) +
    geom_line(colour = "black", lty=1, size = 0.5) +
    geom_line(aes(age, qx.young), colour = "green", lty = 4) +
    geom_line(aes(age, qx.risk), colour = "red", lty = 3, size = 0.5) +
    geom_line(aes(age, qx.adult), colour = "deepskyblue3", lty = 2) +
    geom_point(data = lifeHP, aes(age, qx), shape=1, colour="grey50") +
    ylim(0,0.5) +
    ylab(expression("Mortality (q" [x]* ")")) + xlab("Age") +
    ggtitle("") +
    theme(panel.background =  element_rect(fill = NA, colour = "black", size = 0.5), 
          legend.title = element_blank(), legend.position = "none")
print(HPcurves)


#### LESLIE MATRICES ####

# Set maturity vector: Maturity at age 9 with pregnancy rate of 40% 
mat <- c(0,0,0,0,0,0,0,0,0, rep(1/2.5, 27))


## Leslie Matrix with the Total mortality estimated with the Heligman-Pollard model

# Get the median, low and high 90% credible intervals of the Total Heligman-Pollar mortality
TotalMs <- HP.CI(HPout = modHP, age = 0:35, CI = 90, M = "total")

# Contruct life tables for the median and credible intervals 
TotHP.life <- Est.life.tab(Est.qx = TotalMs$Med, age = 0:35, n = 1000)
TotHP.life.L <- Est.life.tab(Est.qx = TotalMs$Mlo, age = 0:35, n = 1000)
TotHP.life.H <- Est.life.tab(Est.qx = TotalMs$Mhi, age = 0:35, n = 1000)

# Contruct Leslie life tables for the median and credible intervals
Tot.life <- with(TotHP.life, life.Leslie(x = 0:35, nKx = nx, nDx = dx,
                                 type = "cohort", iwidth = 1,
                                 width12 = c(1, 1)))
Tot.lifeL <- with(TotHP.life.L, life.Leslie(x = 0:35, nKx = nx, nDx = dx,
                                   type = "cohort", iwidth = 1,
                                   width12 = c(1,1)))
Tot.lifeH <- with(TotHP.life.H, life.Leslie(x = 0:35, nKx = nx, nDx = dx,
                                   type="cohort", iwidth = 1,
                                   width12  =c(1,1)))

# Contruct Leslie matrices for the median and credible intervals
Tot.A <- Leslie.matrix(lx = Tot.life$nLx, mx = mat, infant.class = FALSE,
            one.sex = TRUE, SRB = 1, L = FALSE, peryear = 1)
Tot.AL <- Leslie.matrix(lx = Tot.lifeL$nLx, mx = mat, infant.class = FALSE,
             one.sex = TRUE, SRB = 1, L = FALSE, peryear = 1)
Tot.AH <- Leslie.matrix(lx = Tot.lifeH$nLx, mx = mat, infant.class = FALSE,
             one.sex = TRUE, SRB = 1, L = FALSE, peryear = 1)

# Predict Leslie projections 
N.tot <- Leslie.pred(A = Tot.A, no = TotHP.life$nx, tmax = 100, pop.sum = TRUE)
NL.tot <- Leslie.pred(A = Tot.AL, no = TotHP.life.L$nx, tmax = 100, pop.sum = TRUE)
NH.tot <- Leslie.pred(A = Tot.AH, no = TotHP.life.H$nx, tmax = 100, pop.sum = TRUE)

# Eigenvalues (medians and credible intervals)
# Population growth
Tot.Aea <- eigen.analysis(Tot.A)
Tot.Aea$lambda1 # 0.9618473
Tot.AeaL <- eigen.analysis(Tot.AL)
Tot.AeaL$lambda1 # 0.9841386
Tot.AeaH <- eigen.analysis(Tot.AH)
Tot.AeaH$lambda1 # 0.9345727
# Net production
Tot.Aea$rho # 1.176908
Tot.AeaL$rho # 1.166733
Tot.AeaH$rho # 1.174052
# Generation time
gen.time(Tot.A, peryear = 1) # 16.44523
gen.time(Tot.AL, peryear = 1) # 17.08953
gen.time(Tot.AH, peryear = 1) # 15.997


## Leslie Matrix with the Natural mortality estimated with the Heligman-Pollard model

# Get the median, low and high 90% credible intervals of the Natural Heligman-Pollar mortality
NaturalMs <- HP.CI(HPout = modHP, age = 0:35, CI = 90, M = "natural")

# Contruct life tables for the median and credible intervals 
NatHP.life <- Est.life.tab(Est.qx = NaturalMs$Med, age = 0:35, n = 1000)
NatHP.life.L <- Est.life.tab(Est.qx = NaturalMs$Mlo, age = 0:35, n = 1000)
NatHP.life.H <- Est.life.tab(Est.qx = NaturalMs$Mhi, age = 0:35, n = 1000)

# Contruct Leslie life tables for the median and credible intervals
Nat.life <- with(NatHP.life, life.Leslie(x = age, nKx = nx, nDx = dx,
                                    type = "cohort", iwidth = 1,
                                    width12 = c(1,1)))
Nat.lifeL <- with(NatHP.life.L, life.Leslie(x = age, nKx = nx, nDx = dx,
                                      type = "cohort", iwidth = 1,
                                      width12 = c(1,1)))
Nat.lifeH <- with(NatHP.life.H, life.Leslie(x = age, nKx = nx, nDx = dx,
                                      type = "cohort", iwidth = 1,
                                      width12 = c(1,1)))

# Contruct Leslie matrices for the median and credible intervals
Nat.A <- Leslie.matrix(lx = Nat.life$nLx, mx = mat, infant.class = FALSE,
            one.sex = TRUE, SRB = 1, L = FALSE, peryear = 1)
Nat.AL <- Leslie.matrix(lx = Nat.lifeL$nLx, mx = mat, infant.class = FALSE,
             one.sex = TRUE, SRB = 1, L = FALSE, peryear = 1)
Nat.AH <- Leslie.matrix(lx = Nat.lifeH$nLx, mx = mat, infant.class = FALSE,
             one.sex = TRUE, SRB = 1, L = FALSE, peryear = 1)

# Predict Leslie projections 
N.nat <- Leslie.pred(A = Nat.A, no = NatHP.life$nx, tmax = 100, pop.sum = TRUE)
NL.nat <- Leslie.pred(A = Nat.AL, no = NatHP.life.L$nx, tmax = 100, pop.sum = TRUE)
NH.nat <- Leslie.pred(A = Nat.AH, no = NatHP.life.H$nx, tmax = 100, pop.sum = TRUE)

# Eigenvalues (medians and credible intervals)
# Population growth
Nat.Aea <- eigen.analysis(Nat.A)
Nat.Aea$lambda1 # 1.030855
Nat.AeaL <- eigen.analysis(Nat.AL)
Nat.AeaL$lambda1 # 1.044063
Nat.AeaH <- eigen.analysis(Nat.AH)
Nat.AeaH$lambda1 # 1.007133
# Net production
Nat.Aea$rho # 1.224137
Nat.AeaL$rho # 1.206632
Nat.AeaH$rho # 1.219181
# Generation time
gen.time(Nat.A, peryear=1) # 16.30796
gen.time(Nat.AL, peryear=1) # 16.92918
gen.time(Nat.AH, peryear=1) # 15.91013

