library(rethinking)

data(Howell1)
d2 <- Howell1[Howell1$age >= 18, ]

dat <- list(N = nrow(d2), height = d2$height)


stan4_1 <- stan(file = "m4_1.stan", data = dat, chains = 2,
                warmup = 1000, iter = 9000, cores = 2)

stan4_1
plot(stan4_1)
traceplot(stan4_1)

samples <- extract(stan4_1)

dens(samples$mu, norm.comp = TRUE)
dens(samples$sigma, norm.comp = TRUE)
hist(samples$mu, breaks = "fd", probability = TRUE, col = 'skyblue')
curve(dnorm(x, mean(samples$mu), sd(samples$mu)), add = TRUE,
      col = "red", lwd =2)


library(rstanarm)

rs4_1 <- stan_glm(height ~ 1,
                  family = gaussian(), 
                  data = d2, 
                  prior_intercept = normal(178, 20, autoscale = FALSE),
                  prior_aux = exponential(1, autoscale = FALSE))

summary(rs4_1, digits = 2)
plot(rs4_1)

posterior_vs_prior(rs4_1)

