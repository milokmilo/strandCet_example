# ################################################################################################ #
#      strandCet: R package for estimating natural and anthropogenic mortality-at-age 
#                        of cetaceans from age-structured strandings.
#                                   Camilo Saavedra
#               Instituto Español de Oceanografía. Centro Oceanográfico de Vigo. 
#                       Subida a Radio Faro Nº50, 36390 Vigo, Spain.
#                         E-mail: camilo.saavedra.penas@gmail.com
# ################################################################################################ #

#### Example ####
  
# In this document can be found the code used to run the examples of the paper 
# titled "*strandCet*: R package for estimating natural and anthropogenic mortality-at-age 
# of cetaceans from age-structured strandings". It should be noted that model fits and 
# estimated parameters could vary slightly from one run to another. General stepwise instructions 
# for the different functions and their options, as well as other examples can be found in 
# the *strandCet* manual <https://cran.r-project.org/web/packages/strandCet/strandCet.pdf>. 
# For further details see the package help files > help(package = "strandCet").
  

### Load necesary packages

# Note that these packages must be previously installed.
library(strandCet)
library(ggplot2)
 

####  Original data ####

# From Stolen and Barlow (2003)
N <- c(36,14,32,15,13,8,4,5,3,3,2,3,3,6,3,8,5,5,8,3,11,4,1,4,3,3,6,2,1,2,0,1,0,1,0,2)
age <- 0:35
data.A <- data.frame(age, N) 
 
# Life table from observed data (original data without bycatch)
life.A <- life.tab(data.A)
print(life.A)


#### Data with the first age class underrepresented ####

# Stolen and Barlow (2003) considered that the number of stranded cetaceans are not negatively 
# biased in the first age classes

# Only one third of the strandings of calves (first age class) were considered in this dataset 
# to simulate a scenario with underrepresentation of calves
Nyoung <- c(N[c(1)]/3, N[-c(1)])
data.B <- data.frame(age, Nyoung) 
 
# Life table from observed data (data with the first age class underrepresented)
life.B <- life.tab(data.B)
print(life.B)


#### Data with theoretical bycaught dolphins added ####

# Stolen and Barlow (2003) considered that the stranded cetaceans belonged to a population 
# not affected by direct anthropogenic mortality

# Bycaught dolphins were added following a gaussian distribution 
# to simulate a population affected by incidental catch
set.seed(123)
Nbyc <- N + c(sort(round(rnorm(6, mean = 10, sd = 5))),
              rev(sort(round(rnorm(6, mean = 10, sd = 5)))), rep(0,24))
data.C <- data.frame(age, Nbyc) 

# Life table from observed data (data with theoretical bycaught dolphins)
life.C <- life.tab(data.C)
print(life.C, row.names = FALSE)
 

#### Siler Model (original data without bycatch) ####

# In the simplest scenario, there are no negative biases in the number of stranded animals 
# of the first age classes and we are interested in calculating only the total mortality 
# (since it is considered that the population is not affected by major anthropogenic pressures)

# Siler model (original data without bycatch)
modSI <-  Si.mod(data.A)
 
# Predict the Siler model (original data without bycatch)
predSI <- Si.pred(data.A, modSI)
 
# Life table of the Siler model (original data without bycatch)
life.Siler <- Est.life.tab(Est.qx = predSI$qx.tot, age = 0:35, n = 1000)
print(life.Siler)
 
# Siler plot (original data without bycatch)
Sicurves <- ggplot(predSI, aes(age, qx.tot)) +
  geom_line(colour = "red", lty = 1, size = 0.5) +
  geom_point(data = life.A, aes(age, qx), shape = 1, colour = "grey50") +
  ylim(0, 0.5) +
  ylab(expression("Mortality (q" [x]* ")")) + xlab("Age") +
  ggtitle("") +
  theme(panel.background =  element_rect(fill = NA, colour = "black", size = 0.5), 
        legend.title = element_blank(), legend.position = "none")
print(Sicurves)
 

#### Siler Model (data with the first age class underrepresented) ####

# The first age classes are usually underrepresented in strandings. 
# One option is to remove the ages considered negatively biased fitting the model 
# to the remaining ages. If this bias/underrepresentation is ignored, we can get wrong fits, 
# even with lower mortality rates in the first age classes than in the juvenile ages.

# Siler model (data with the first age class underrepresented)
modSI.R <-  Si.mod(data.B, rm = 1)
modSI <-  Si.mod(data.B)
 
# Predict the Siler model (data with the first age class underrepresented)
predSI.R <- Si.pred(data.B, modSI.R, rm = 1)
predSI <- Si.pred(data.B, modSI)
 
# Life table of the Siler model (data with the first age class underrepresented)
life.Siler <- Est.life.tab(Est.qx = predSI.R$qx.tot, age = 0:35, n = 1000)
print(life.Siler)
 
# Siler plot (data with the first age class underrepresented)
Sicurves <- ggplot(predSI.R, aes(age, qx.tot)) +
  geom_line(colour = "red", lty = 1, size = 0.5) +
  geom_line(data = predSI, colour = "black", lty = 2, size = 0.5) +
  geom_point(data = life.B, aes(age, qx), shape = 1, colour = "grey50") +
  geom_point(data = life.B[c(1),], aes(age, qx), shape = 16, colour="grey50") +
  ylim(0, 0.5) +
  ylab(expression("Mortality (q" [x]* ")")) + xlab("Age") +
  ggtitle("") +
  theme(panel.background =  element_rect(fill = NA, colour = "black", size = 0.5), 
        legend.title = element_blank(), legend.position = "none")
print(Sicurves)
 

#### Heligman-Pollard Model (data with theoretical bycaught dolphins) ####

# Mortality-at-age of a population affected by high levels of bycatch 
# (which mainly affects young ages) should not be fitted with a Siler model. 
# Mortality rates may not show a tipical U-shape, so the fit could not be adequate. 
# In addition, it is only posible to calculate the total mortality. In these cases, 
# the adapted Heligman-Pollard model is more appropriate.

# Siler model (data with theoretical bycaught dolphins)
modSI <-  Si.mod(data.C)

# Predict the Siler model (data with theoretical bycaught dolphins)
predSI <- Si.pred(data.C, modSI)
 
# Select non-informative priors for the HP model (data with theoretical bycaught dolphins)
priors <- data.frame(priors.lo = c(0,0,0,0,0,0,1,0,1),
                     priors.hi = c(1,10,1,0.01,0.5,10,15,0.01,1.5))

# Compile priors in the required format (data with theoretical bycaught dolphins)
q0 <- HP.priors(pri.lo = priors$priors.lo,
                pri.hi = priors$priors.hi,
                theta.dim = 9)
 
# Run the Heligman-Pollard model (data with theoretical bycaught dolphins)
modHP <- HP.mod(prior = q0, lifeTab = life.C,
                K = 10, d = 10, B = 500, CI = 90)
print(modHP$out, row.names = c("A", "B", "C", "I", "D", "E", "F", "G", "H"))
 
# Predict the Heligman-Pollard model (data with theoretical bycaught dolphins)
predHP <- HP.pred(life = life.C, HPout = modHP, age = age)
print(predHP, row.names = FALSE)
 
# Heligman-Pollard plot (data with theoretical bycaught dolphins)
HPcurves <- ggplot(predHP, aes(age, qx.tot)) +
  geom_line(colour = "black", lty=1, size = 0.8) +
  geom_line(data = predSI, colour = "black", lty=1, size = 0.5) +
  geom_line(aes(age, qx.young), colour = "green", lty = 4) +
  geom_line(aes(age, qx.risk), colour = "red", lty = 3, size = 0.5) +
  geom_line(aes(age, qx.adult), colour = "deepskyblue3", lty = 2) +
  geom_point(data = life.C, aes(age, qx), shape=1, colour="grey50") +
  ylim(0,0.65) +
  ylab(expression("Mortality (q" [x]* ")")) + xlab("Age") +
  ggtitle("") +
  theme(panel.background =  element_rect(fill = NA, colour = "black", size = 0.5), 
        legend.title = element_blank(), legend.position = "none")
#png("aHPmodel.png", height = 10, width = 20, units = 'cm', res = 900)
print(HPcurves)
#dev.off()


#### Leslie Matrices (data with theoretical bycaught dolphins) ####

# Three generation time was used for projections following IUCN recommendations.
# Since the generation time of the bottlenose dolphins (*Tursiops truncatus*) used in this analysis 
# is about 21 years (Taylor et al. 2007) 65 years were used as an approximation of three generations

# Set a maturity vector:
# Maturity at age 7 (first parturition at age 8) with pregnancy rate of 40% (Wells et al. 1987)
mat <- c(0,0,0,0,0,0,0,0, rep(0.4, 28))


#### Leslie Matrix using the total mortality estimated with the Heligman-Pollard model ####

# Get the median, low and high 90% credible intervals 
# of the Natural Heligman-Pollar mortality (using the total mortality)
TotalMs <- HP.CI(HPout = modHP, age = 0:35, CI = 90, M = "total")
 
# Contruct life tables for the median and credible intervals (using the total mortality)
TotHP.life <- Est.life.tab(Est.qx = TotalMs$Med, age = 0:35, n = 1000)
TotHP.life.L <- Est.life.tab(Est.qx = TotalMs$Mlo, age = 0:35, n = 1000)
TotHP.life.H <- Est.life.tab(Est.qx = TotalMs$Mhi, age = 0:35, n = 1000)
 
# Contruct Leslie life tables for the median and credible intervals (using the total mortality)
Tot.life <- with(TotHP.life, life.Leslie(x = 0:35, nKx = nx, nDx = dx,
                                         type = "cohort", iwidth = 1,
                                         width12 = c(1, 1)))
Tot.lifeL <- with(TotHP.life.L, life.Leslie(x = 0:35, nKx = nx, nDx = dx,
                                            type = "cohort", iwidth = 1,
                                            width12 = c(1,1)))
Tot.lifeH <- with(TotHP.life.H, life.Leslie(x = 0:35, nKx = nx, nDx = dx,
                                            type="cohort", iwidth = 1,
                                            width12  =c(1,1)))
 
# Contruct Leslie matrices for the median and credible intervals (using the total mortality)
Tot.A <- Leslie.matrix(lx = Tot.life$nLx, mx = mat, infant.class = FALSE,
                       one.sex = TRUE, SRB = 1, L = FALSE, peryear = 1)
Tot.AL <- Leslie.matrix(lx = Tot.lifeL$nLx, mx = mat, infant.class = FALSE,
                        one.sex = TRUE, SRB = 1, L = FALSE, peryear = 1)
Tot.AH <- Leslie.matrix(lx = Tot.lifeH$nLx, mx = mat, infant.class = FALSE,
                        one.sex = TRUE, SRB = 1, L = FALSE, peryear = 1)
 
# Predict Leslie matrices (using the total mortality)
N.tot <- Leslie.pred(A = Tot.A, no = TotHP.life$nx, tmax = 65, pop.sum = TRUE)
NL.tot <- Leslie.pred(A = Tot.AL, no = TotHP.life.L$nx, tmax = 65, pop.sum = TRUE)
NH.tot <- Leslie.pred(A = Tot.AH, no = TotHP.life.H$nx, tmax = 65, pop.sum = TRUE)
 
# Population growth
Tot.Aea <- eigen.analysis(Tot.A); Tot.Aea$lambda1
Tot.AeaL <- eigen.analysis(Tot.AL); Tot.AeaL$lambda1
Tot.AeaH <- eigen.analysis(Tot.AH); Tot.AeaH$lambda1
 
# Net production
Tot.Aea$rho
Tot.AeaL$rho
Tot.AeaH$rho
 
# Generation time
gen.time(Tot.A, peryear = 1)
gen.time(Tot.AL, peryear = 1)
gen.time(Tot.AH, peryear = 1)
 

#### Leslie Matrix using the natural mortality estimated with the Heligman-Pollard model ####

# Get the median, low and high 90% credible intervals 
# of the Natural Heligman-Pollar mortality (using the natural mortality)
NaturalMs <- HP.CI(HPout = modHP, age = 0:35, CI = 90, M = "natural")
 
# Contruct life tables for the median and credible intervals (using the natural mortality)
NatHP.life <- Est.life.tab(Est.qx = NaturalMs$Med, age = 0:35, n = 1000)
NatHP.life.L <- Est.life.tab(Est.qx = NaturalMs$Mlo, age = 0:35, n = 1000)
NatHP.life.H <- Est.life.tab(Est.qx = NaturalMs$Mhi, age = 0:35, n = 1000)
 
# Contruct Leslie life tables for the median and credible intervals (using the natural mortality)
Nat.life <- with(NatHP.life, life.Leslie(x = age, nKx = nx, nDx = dx,
                                         type = "cohort", iwidth = 1,
                                         width12 = c(1,1)))
Nat.lifeL <- with(NatHP.life.L, life.Leslie(x = age, nKx = nx, nDx = dx,
                                            type = "cohort", iwidth = 1,
                                            width12 = c(1,1)))
Nat.lifeH <- with(NatHP.life.H, life.Leslie(x = age, nKx = nx, nDx = dx,
                                            type = "cohort", iwidth = 1,
                                            width12 = c(1,1)))

# Contruct Leslie matrices for the median and credible intervals (using the natural mortality)
Nat.A <- Leslie.matrix(lx = Nat.life$nLx, mx = mat, infant.class = FALSE,
                       one.sex = TRUE, SRB = 1, L = FALSE, peryear = 1)
Nat.AL <- Leslie.matrix(lx = Nat.lifeL$nLx, mx = mat, infant.class = FALSE,
                        one.sex = TRUE, SRB = 1, L = FALSE, peryear = 1)
Nat.AH <- Leslie.matrix(lx = Nat.lifeH$nLx, mx = mat, infant.class = FALSE,
                        one.sex = TRUE, SRB = 1, L = FALSE, peryear = 1)

# Predict Leslie matrices (using the natural mortality)
N.nat <- Leslie.pred(A = Nat.A, no = NatHP.life$nx, tmax = 65, pop.sum = TRUE)
NL.nat <- Leslie.pred(A = Nat.AL, no = NatHP.life.L$nx, tmax = 65, pop.sum = TRUE)
NH.nat <- Leslie.pred(A = Nat.AH, no = NatHP.life.H$nx, tmax = 65, pop.sum = TRUE)
 
# Population growth
Nat.Aea <- eigen.analysis(Nat.A); Nat.Aea$lambda1
Nat.AeaL <- eigen.analysis(Nat.AL); Nat.AeaL$lambda1
Nat.AeaH <- eigen.analysis(Nat.AH); Nat.AeaH$lambda1
 
# Net production
Nat.Aea$rho
Nat.AeaL$rho
Nat.AeaH$rho

# Generation time
gen.time(Nat.A, peryear = 1)
gen.time(Nat.AL, peryear = 1)
gen.time(Nat.AH, peryear = 1)


#### Bibliography ####

# Stolen M., Barlow J. 2003. A model life table for bottlenose dolphins (*Tursiops truncatus*) 
# from the Indian River Lagoon System, Florida, USA. *Marine mammal science* 19:630–649. 

# Taylor, B.L., Chivers, S.J., Larese, J., Perrin, W.F., 2007. Generation length and percent 
# mature estimates for IUCN assessments of cetaceans. NOAA Technical Memorandum SFSC.

# Wells, RS., Scott, MD., & Irvine, AB. (1987). The social structure of free-ranging bottlenose 
# dolphins. In: HH. Genoways (Ed.), Current Mammalogy (pp. 247–305). New York: Plenum Press.

