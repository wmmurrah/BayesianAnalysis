#----------EXAMPLE 6.3: MODEL COMPARISON----------------#

# Model Comparison: Background variables only

BGModel_inf <- MCMCregress(rcomb1~gender+native+ slang+ ESCS
       ,data=datafile9,marginal.likelihood="Chib95",mcmc=10000,
	b0=c(470.9, 26.3, 4.7, 23.3, 39.9 ),
	B0=c( 0.0185, 0.0952, 0.0151, 0.0222, 0.3541 ))
plot(BGModel_inf)
dev.off()
summary(BGModel_inf)

# Model Comparison: Attitudinal variables only

ATTModel_inf <- MCMCregress(rcomb1~JOYREAD+ DIVREAD,
	data=datafile9,marginal.likelihood="Chib95",mcmc=10000,
	b0=c( 505.4, 27.2, 8.4 ),B0=c(0.3643, 0.3147, 0.3497))
plot(ATTModel_inf)
dev.off()
summary(ATTModel_inf)

# Model Comparison: Learning strategies variables only

LSModel_inf <- MCMCregress(rcomb1~ MEMOR+ ELAB+ CSTRAT,
	data=datafile9,marginal.likelihood="Chib95",mcmc=10000,
	b0=c( 509.7, -24.2, -9.8, 38.9),B0=c( 0.3327, 0.1829, 0.1848, 0.1563))
plot(LSModel_inf)
dev.off()
summary(LSModel_inf)

# Calculation of Bayes Factors

bf <- BayesFactor(BGModel_inf, ATTModel_inf, LSModel_inf, FullModel_inf)
print(bf)

