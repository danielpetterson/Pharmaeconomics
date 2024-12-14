# rm(list=ls())
# setwd("C:\\Users\\Salmasi\\Dropbox\\Cattolica\\Didattica\\Pharmaeconomics and HTA\\Case studies\\Lab sessions\\R-scripts")

# Loads the package to run OpenBUGS or JAGS from R
library(R2OpenBUGS)
library(R2jags)

# Launches the file Utils.R which contains useful functions used throughout this script
source("Utils.R")
# Loads the data into R (assumes the file is stired in the working directory - if not the full path can be provided)  
source("LoadData_exams.R")
# Sets the working directory
# setwd("C:\\Users\\Salmasi\\Dropbox\\Cattolica\\Didattica\\Pharmaeconomics and HTA\\Case studies\\Lab sessions\\R-scripts")
working.dir <- getwd()
# setwd("C:/Users/Salmasi/Dropbox/Cattolica/Didattica/Pharmaeconomics and HTA")

# Defines the data list to be passed to BUGS/JAGS
data <- list("N","a.phi","b.phi","mu.rho","tau.rho","a.beta","b.beta","a.gamma","b.gamma","mu.omega","tau.omega","mu.psi","tau.psi","N.outcomes","N.resources","mu.lambda","tau.lambda","a.xi","b.xi","a.eta","b.eta","a.delta")

# Defines the file with the model code
filein <- "vaccine_mod.R"

# Defines the quantities to be monitored (stored)
params <- c("beta","phi","omega","rho","Infected","GP","Repeat.GP","Pneumonia","Hospital","Death","Mild.Compl","Trt","Adverse.events","n","gamma","delta","psi","lambda","pi","xi","eta")

# Generates the initial values
inits <- function(){
  list(phi=runif(1),beta=runif(N.outcomes,0,1),rho=c(NA,runif(1)),gamma=runif(2,0,1),delta=rpois(1,2),omega=c(runif(1),NA,NA,runif(1),NA,runif(2,0,1)),psi=runif(N.resources,0,10),lambda=runif(1),eta=runif(1),xi=runif(1))
}

# Defines the number of iteration, burn-in and thinning, and runs BUGS or JAGS

n.iter <- 100000
n.burnin <- 9500
n.thin <- floor((n.iter-n.burnin)/500)

# 1. This runs OpenBUGS - Very ineffiecient, takes 4 days to perform simulations
# vaccine_bugs <- bugs(data, inits, params, model.file = filein, n.chains = 2, n.iter, n.burnin, n.thin, DIC = FALSE, working.directory = working.dir, debug = FALSE)

# 2. This runs JUGS
vaccine_jags <- jags(data, inits, params, model.file = filein, n.chains = 2, n.iter, n.burnin, n.thin, DIC = FALSE, working.directory = working.dir, progress.bar = "text")
#vaccine_jags <- jags(data, inits = NULL, params, model.file = filein, n.chains = 2, n.iter, n.burnin, n.thin, DIC = FALSE, working.directory = working.dir, progress.bar = "text")

# Prints the summary stats and attaches the results to the R workspace
print(vaccine_jags, digits = 3, intervals=c(0.025,0.975))

# In OpenBUGS
#attach.bugs(vaccine)
# In JAGS
attach.jags(vaccine_jags)

# traceplot
# traceplot(vaccine_jags)

save.image("RunMCMC.RData")

#load("RunMCMC.RData")

# Economic model

# load("http://www.statistica.it/gianluca/BCEABook/vaccine.RData")

# Compute effectiveness in QALYs lost for both strategies

# Infected[h,t,c] is an array of h matrices matrices of dimension 2X2

Infected[1,,] # time periods, t, are the rows of each matrix
              # non vaccinated/vaccinated subjects, v, are the columns of each matrix

QALYs.inf <- QALYs.pne <- QALYs.hosp <- QALYs.adv <- QALYs.death <- matrix(0,n.sims,2)

for (t in 1:2) {
  QALYs.inf[,t] <- ((Infected[,t,1] + Infected[,t,2])*omega[,1]/350)
  QALYs.pne[,t] <- ((Pneumonia[,t,1] + Pneumonia[,t,2])*omega[,4]/350)
  QALYs.hosp[,t] <- ((Hospital[,t,1] + Hospital[,t,2])*omega[,5]/350)
  QALYs.death[,t] <- ((Death[,t,1] + Death[,t,2])*omega[,6]/350)
  QALYs.adv[,t] <- (Adverse.events[,t,2]*omega[,7]/350)
}

# population average measure of effectiveness
e <- -(QALYs.inf + QALYs.pne + QALYs.adv + QALYs.hosp + QALYs.death)/N

summary(e)

# Compute costs for both strategies

cost.GP <- cost.hosp <- cost.vac <- cost.time.vac <- cost.time.off <- cost.trt1 <- cost.trt2 <- cost.otc <- cost.travel <- matrix(0, n.sims, 2)

for (t in 1:2) {
  cost.GP[,t] <- (GP[,t,1] + GP[,t,2] + Repeat.GP[,t,1] + Repeat.GP[,t,2])*psi[,1]
  cost.hosp[,t] <- (Hospital[,t,1] + Hospital[,t,2])*psi[,2]
  cost.vac[,t] <- n[,t,2]*psi[,3]
  cost.time.vac[,t] <- n[,t,2]*psi[,4]
  cost.time.off[,t] <- (Infected[,t,1] + Infected[,t,2])*psi[,5]*eta*lambda #psi[,5] = Cost of days work absence due to influenza, eta = Being off work, lanbda = lenght of absence from work
  cost.trt1[,t] <- (GP[,t,1] + GP[,t,2])*gamma[,1]*psi[,6]*delta # gamma = prob to take AVD, psi[,6] = cost, delta = number 
  cost.trt2[,t] <- (Repeat.GP[,t,1] + Repeat.GP[,t,2])*gamma[,2]*psi[,6]*delta
  cost.otc[,t] <- (Infected[,t,1] + Infected[,t,2])*psi[,7]*xi # psi[,7] = cost of OTC, xi = prob to take OTC
  cost.travel[,t] <- n[,t,2]*psi[,8]
}

# average population cost
c <- (cost.GP + cost.hosp + cost.vac + cost.time.vac + cost.time.off + cost.trt1 + cost.trt2 + cost.travel + cost.otc)/N

summary(c)
summary(e)

# Cost-Effectiveness analysis

# install.packages("BCEA")
library(BCEA)
treats <- c("Status quo", "Vaccination")
m <- bcea(e,c,ref=2, plot = TRUE)

summary(m)
summary(m, wtp=10000)
summary(m, wtp=30000)

jpeg(file="ceplane1.jpg", width=700, height=800)
ceplane.plot(m, xlab = "Difference in QALYs", ylab = "Difference in costs (Pounds)", title = "C/E plane", cex = 2.5)
dev.off()
jpeg(file="ceplane2.jpg", width=700, height=800)
ceplane.plot(m, wtp = 30000, xlab = "Difference in QALYs", ylab = "Difference in costs (Pounds)", title = "C/E plane")
dev.off()
jpeg(file="ceplane3.jpg", width=700, height=800)
ceplane.plot(m, wtp = 10000, xlab = "Difference in QALYs", ylab = "Difference in costs (Pounds)", title = "C/E plane")
dev.off()
jpeg(file="ibplot_30.jpg", width=700, height=800)
ib.plot(m, wtp=30000)
dev.off()
jpeg(file="ibplot_50.jpg", width=700, height=800)
ib.plot(m, wtp=50000)
dev.off()
jpeg(file="ibplot_10.jpg", width=700, height=800)
ib.plot(m, wtp=10000)
dev.off()
jpeg(file="eibplot.jpg", width=700, height=800)
eib.plot(m)
dev.off()
jpeg(file="contour_1.jpg", width=700, height=800)
contour(m)
dev.off()
jpeg(file="contour_2.jpg", width=700, height=800)
contour2(m, wtp=30000)
dev.off()

# Probabilistic Sensitivity Ananlysis (PSA)

plot(m,
     graph="ggplot2",                                      # use ggplot2
     theme=theme(plot.title=element_text(size=rel(1.25))), # theme elements must have a name
     ICER.size=1.5,                                        # hidden option in ceplane.plot
     size=rel(2.5)                                         # modifies the size of k= labels
)                                                       #  in ceplane.plot and eib.plot

table <- sim_table(m, wtp = 30000)
head(table$Table)

jpeg(file="ceplane_ceac.jpg", width=700, height=800)
par(mfrow=c(2,2))
ceplane.plot(m, wtp = 1000, xlab = "Difference in QALYs", ylab = "Difference in costs (Pounds)", title = "C/E plane")
ceplane.plot(m, wtp = 10000, xlab = "Difference in QALYs", ylab = "Difference in costs (Pounds)", title = "C/E plane")
ceplane.plot(m, wtp = 30000, xlab = "Difference in QALYs", ylab = "Difference in costs (Pounds)", title = "C/E plane")
ceac.plot(m)
dev.off()

jpeg(file="ceac.jpg", width=700, height=800)
ceac.plot(m)
dev.off()

jpeg(file="evi.jpg", width=700, height=800)
evi.plot(m)
dev.off()

inp <- createInputs(vaccine_jags)
names(inp)

EVPPI <- evppi(m,c("beta.1.", "beta.2."),inp$mat)

jpeg(file="evppi.jpg", width=700, height=800)
plot(EVPPI)
dev.off()

EVPPI$evppi[which(EVPPI$k==30000)]

EVPPI <- evppi(m,c("beta.1.", "beta.2."),inp$mat,residuals=T)

jpeg(file="evppi_res.jpg", width=700, height=800)
diag.evppi(EVPPI,m,plot_type ="residuals")
dev.off()

jpeg(file="evppi_qqplot.jpg", width=700, height=800)
diag.evppi(EVPPI,m,plot_type="qqplot")
dev.off()

EVPPI.so <- evppi(m,c("beta.1.", "beta.2."),inp$mat,method="so",n.blocks=50)
EVPPI.sad <- evppi(m,c("beta.1.", "beta.2."),inp$mat,method="sad",n.seps=1)

jpeg(file="evppi_so.jpg", width=700, height=800)
plot(EVPPI.so)
dev.off()

jpeg(file="evppi_sad.jpg", width=700, height=800)
plot(EVPPI.sad)
dev.off()

jpeg(file="info_rank.jpg", width=700, height=800)
info.rank(m, inp, graph = "base")
dev.off()
