install.packages(c("R2jags",
                   "BCEA",
                   "splancs",
                   "ggplot2",
                   "heemod",
                   "diagram",
                   "R2OpenBUGS",
                   "R2jags"))
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"),
                 dep=TRUE)

# Change wd to local Lab Session path
wd <- '/Users/danielpetterson/GitHub/Pharmaeconomics/Lab Sessions'
setwd(wd)

source("Markov_model.R")
source("Statistical_and_Economic_model_Vaccine.R")
