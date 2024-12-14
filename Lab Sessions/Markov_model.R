#https://cran.r-project.org/web/packages/heemod/vignettes/c_homogeneous.html

# rm(list=ls())
# setwd("C:\\Users\\Salmasi\\Dropbox\\Cattolica\\Didattica\\Pharmaeconomics and HTA\\Case studies\\Lab sessions\\R-scripts")

library(heemod)
library(ggplot2)

param <- define_parameters(
  rr = .509,
  dr = .06,
  
  p_AA_mono = .721,
  p_AB_mono = .202,
  p_AC_mono = .067,
  p_AD_mono = .010,
  
  p_BC_mono = .407,
  p_BD_mono = .012,
  
  p_CD_mono = .250,
  
  
  p_AB_comb = p_AB_mono*rr,
  p_AC_comb = p_AC_mono*rr,
  p_AD_comb = p_AD_mono*rr,
  
  p_BC_comb = p_BC_mono*rr,
  p_BD_comb = p_BD_mono*rr,
  
  p_CD_comb = p_CD_mono*rr,
  
  p_AA_comb = 1 - (p_AB_comb + p_AC_comb + p_AD_comb),
  
  cost_zido = 2278,
  cost_lami = 2086,
  
  cost_A = 2756,
  cost_B = 3052,
  cost_C = 9007
)

mat_trans_mono <- define_transition(
  C, p_AB_mono, p_AC_mono, p_AD_mono,
  0,         C,         p_BC_mono, p_BD_mono,
  0,         0,         C,         p_CD_mono,
  0,         0,         0,         1
)

mat_trans_comb <- define_transition(
  C, p_AB_comb, p_AC_comb, p_AD_comb,
  0,         C,         p_BC_comb, p_BD_comb,
  0,         0,         C,         p_CD_comb,
  0,         0,         0,         1
)

# Just for the plot

mat_trans_mono_plot <- define_transition(
  .722, .202, .067, .010,
  0,        0.581,   .407, 0.012,
  0,         0,         .75,         .25,
  0,         0,         0,         1
)

plot(mat_trans_mono_plot)

plot(mat_trans_mono)
plot(mat_trans_comb)

state_A <- define_state(
  cost_health = cost_A,
  cost_drugs = dispatch_strategy(
    mono = cost_zido,
    comb = cost_zido + cost_lami
  ),
  cost_total = discount(cost_health + cost_drugs, dr),
  life_year = 1
)
state_B <- define_state(
  cost_health = cost_B,
  cost_drugs = dispatch_strategy(
    mono = cost_zido,
    comb = cost_zido + cost_lami
  ),
  cost_total = discount(cost_health + cost_drugs, dr),
  life_year = 1
)
state_C <- define_state(
  cost_health = cost_C,
  cost_drugs = dispatch_strategy(
    mono = cost_zido,
    comb = cost_zido + cost_lami
  ),
  cost_total = discount(cost_health + cost_drugs, dr),
  life_year = 1
)
state_D <- define_state(
  cost_health = 0,
  cost_drugs = 0,
  cost_total = discount(cost_health + cost_drugs, dr),
  life_year = 0
)

strat_mono <- define_strategy(
  transition = mat_trans_mono,
  state_A,
  state_B,
  state_C,
  state_D
)

strat_comb <- define_strategy(
  transition = mat_trans_comb,
  state_A,
  state_B,
  state_C,
  state_D
)

res_mod <- run_model(
  mono = strat_mono,
  comb = strat_comb,
  parameters = param,
  cycles = 50,
  cost = cost_total,
  effect = life_year
)

summary(res_mod, threshold = c(1000, 5000, 6000, 1e4))

plot(res_mod, type = "counts", panel = "by_strategy") +
  xlab("Time") + theme_bw() + scale_color_brewer(name = "State",palette = "Set1")

plot(res_mod, type = "counts", panel = "by_state") +
  xlab("Time") + theme_bw() + scale_color_brewer(name = "Strategy", palette = "Set1")

plot(res_mod, type = "values", panel = "by_value", free_y = TRUE) + 
  xlab("Time") + theme_bw() + scale_color_brewer(name = "Strategy", palette = "Set1")

# DSA

se <- define_dsa(
  rr, .4, .6,
  
  cost_zido, 1500, 3000,
  cost_lami, 1500, 3000,

  dr, .04, .08
)

res_dsa <- run_dsa(
  model = res_mod,
  dsa = se
)

plot(res_dsa, strategy = "comb", result = "icer", type = "difference")
plot(res_dsa, strategy = "comb", result = "icer", type = "difference", limits_by_bars = FALSE)

# PSA

rsp <- define_psa(
  rr ~ lognormal(mean = .509, sdlog = .173),
  
  cost_A ~ gamma(mean = 2756, sd = sqrt(2756)),
  cost_B ~ gamma(mean = 3052, sd = sqrt(3052)),
  cost_C ~ gamma(mean = 9007, sd = sqrt(9007)),
  
#  p_CD_mono ~ binomial(prob = .25, size = 40), # mean proportion and the size of the sample used to estimate that proportion.
  p_CC_mono + p_CD_mono ~ multinomial(750, 250), # mean # multinomial distributions are declared with the number of individuals in each group in the sample used to estimate the proportions. 
  p_BB_mono + p_BC_mono + p_BD_mono ~ multinomial(581, 407, 12), # multinomial distributions are declared with the number of individuals in each group in the sample used to estimate the proportions. 
  p_AA_mono + p_AB_mono + p_AC_mono + p_AD_mono ~ multinomial(721, 202, 67, 10) # multinomial distributions are declared with the number of individuals in each group in the sample used to estimate the proportions. 
)

pm <- run_psa(
  model = res_mod,
  psa = rsp,
  N = 100
)

summary(pm, threshold = c(1000, 5000, 6000, 1e4))

plot(pm, type = "ce")

plot(pm, type = "ac", max_wtp = 10000, log_scale = FALSE)
plot(pm, type = "evpi", max_wtp = 10000, log_scale = FALSE)

#plot(pm, type = "cov")
#plot(pm, type = "cov", diff = TRUE, threshold = 5000)

#Parallel computing
#Resampling can be significantly sped up by using parallel computing. This can be done in the following way:
  
#Define a cluster with the use_cluster() functions (i.e. use_cluster(4) to use 4 cores).
#Run the analysis as usual.
#To stop using parallel computing use the close_cluster() function.
#Results may vary depending on the machine, but we found speed gains to be quite limited beyond 4 cores.

e_mm <- cbind(pm$psa$.effect[pm$psa$.strategy_names == "mono"], pm$psa$.effect[pm$psa$.strategy_names == "comb"])
c_mm <- cbind(pm$psa$.cost[pm$psa$.strategy_names == "mono"], pm$psa$.cost[pm$psa$.strategy_names == "comb"]) 

library(BCEA)

mm <- bcea(e_mm, c_mm, ref = 2, plot = T)

summary(mm, wtp = 10000)

