model {
# 1. Define the number of people in each group n[t,v], where t=1,2 is status quo (t=1) vs vaccination (t=2) 
  # and v=1,2 is non vaccinated vs vaccinated
# t=1: If the vaccine is not available, no one will use it 
	# number of vaccinated in the population
	V[1] <- 0
	# number of individuals in the two groups
	n[1,1] <- N - V[1] 	# non vaccinated
	n[1,2] <- V[1]			# vaccinated

# t=2: When the vaccine is available, some people will use it but some people won't
	# number of vaccinated in the population
	V[2] ~ dbin(phi,N)  # Likelihood
	# number of individuals in the two groups
	n[2,1] <- N - V[2]	# non vaccinated
	n[2,2] <- V[2] 			# vaccinated

# 2. Vaccination coverage
	phi ~ dbeta(a.phi,b.phi) # Prior

# 3. Probability of experiencing the clinical outcomes (in total, N.outcomes = 7)
	# 1. Influenza
	# 2. GP visits
	# 3. Minor complications (repeat visit)
	# 4. Major complications (pneumonia)
	# 5. Hospitalisations 
	# 6. Death
	# 7. Adverse events due to vaccination
	for (r in 1:4) {
		beta[r] ~ dbeta(a.beta[r],b.beta[r]) # All priors
	}
	for (r in 5:6) {
		beta[r] ~ dlnorm(a.beta[r],b.beta[r])
	}
	beta[N.outcomes] ~ dbeta(a.beta[N.outcomes],b.beta[N.outcomes])
#	beta[7] ~ dbeta(a.beta[7],b.beta[7]) # adverse events due to vaccination
	
# 4. Vaccine effectiveness in reducing influenza (for v=1, it is obviously 0)
	rho[1] <- 0 # for non vaccinated effectiveness = 0
	rho[2] ~ dlnorm(mu.rho, tau.rho) # Prior - for vaccinated effectiveness > 0

# 5. Probability of influenza infection
	for (t in 1:2) {
		for (v in 1:2) {
		   pi[t,v] <- beta[1]*(1-rho[v]) # Likelihood
		} 
	}

# 6. Number of patients experiencing the events for both interventions & compliance groups
    
	Infected[1,1] ~ dbin(pi[1,1],n[1,1])      # status quo/non vaccinated
	GP[1,1] ~ dbin(beta[2],Infected[1,1])     # All Likelihoods
	Repeat.GP[1,1] ~ dbin(beta[3],GP[1,1])
	Pneumonia[1,1] ~ dbin(beta[4],GP[1,1])
	Hospital[1,1] ~ dbin(beta[5],GP[1,1])
	Death[1,1] ~ dbin(beta[6],GP[1,1])
	Trt[1,1,1] ~ dbin(gamma[1],GP[1,1])
	Trt[2,1,1] ~ dbin(gamma[2],Mild.Compl[1,1])
	Mild.Compl[1,1] <- Repeat.GP[1,1] + Pneumonia[1,1]
	Adverse.events[1,1] <- 0
	
	Infected[1,2] <- 0                   # status quo/vaccinated all 0s
	GP[1,2] <- 0
	Repeat.GP[1,2] <- 0 
	Pneumonia[1,2] <- 0 
	Hospital[1,2] <- 0 
	Death[1,2] <- 0 
	Trt[1,1,2] <- 0 
	Trt[2,1,2] <- 0 
	Mild.Compl[1,2] <- Repeat.GP[1,2] + Pneumonia[1,2]
	Adverse.events[1,2] <- 0
	
	for (v in 1:2) {
	Infected[2,v] ~ dbin(pi[2,v],n[2,v])        # vaccine available: non vaccinated/vaccinated
	GP[2,v] ~ dbin(beta[2],Infected[2,v])
	Repeat.GP[2,v] ~ dbin(beta[3],GP[2,v])
	Pneumonia[2,v] ~ dbin(beta[4],GP[2,v])
	Hospital[2,v] ~ dbin(beta[5],GP[2,v])
	Death[2,v] ~ dbin(beta[6],GP[2,v])
	Trt[1,2,v] ~ dbin(gamma[1],GP[2,v])
	Trt[2,2,v] ~ dbin(gamma[2],Mild.Compl[2,v])
	Mild.Compl[2,v] <- Repeat.GP[2,v] + Pneumonia[2,v]
	}
	
	Adverse.events[2,1] <- 0
	Adverse.events[2,2]	~ dbin(beta[N.outcomes],n[2,2])
	
# 7. Probability of experiencing other events (impacts on costs and QALYs/QALDs)
# Treatment with antibiotics after GP visit
	
	for (i in 1:2) {
	  gamma[i] ~ dbeta(a.gamma[i],b.gamma[i]) # Prior
	}
	# Number of prescriptions of antivirals
	delta ~ dpois(a.delta) # Prior
	# Taking OTC
	xi ~ dbeta(a.xi,b.xi) # Prior
	# Being off work
	eta ~ dbeta(a.eta,b.eta) # Prior
	# Length of absence from work from influenza
	lambda ~ dlnorm(mu.lambda,tau.lambda) # Prior

# 8. Costs of clinical resources (N.resources = 8)
	# 1. Cost of GP visit
	# 2. Cost of hospital episode
	# 3. Cost of vaccination
	# 4. Cost of time to receive vaccination
	# 5. Cost of days work absence due to influenza
	# 6. Cost of antiviral drug
	# 7. Cost of OTC treatments
	# 8. Cost of travel to receive vaccination
	for (r in 1:N.resources) {
		psi[r] ~ dlnorm(mu.psi[r],tau.psi[r]) # Prior
	}

# 9. Quality of life adjusted days/years loss
	# 1. Influenza infection
	# 2. GP visits (no QALD/QALY loss)
	# 3. Minor complications (repeat visit, no QALD/Y loss)
	# 4. Major complications (pneumonia)
	# 5. Hospitalisations (same QALD/Y loss as pneumonia)
	# 6. Death
	# 7. Adverse events due to vaccination
	omega[1] ~ dlnorm(mu.omega[1],tau.omega[1]) # Prior
	omega[2] <- 0; omega[3] <- 0; # Prior
	for (r in 4:N.outcomes) {
		omega[r] ~ dlnorm(mu.omega[r],tau.omega[r]) # Prior
	} 

}			