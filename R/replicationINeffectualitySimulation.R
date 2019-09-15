
# Parameters
N       = 100                  # Number of labs
b       = 0.1                  # base rate of true hypothesis
r0      = c(0, 0.01, 0.2, 0.5) # initial replication rate of all labs
e0      = 75                   # initial effort of all labs
w0      = 0.80                 # initial power of all labs
eta     = 0.2                  # influence of effort on productivity
Crpos  = 1                    # probability of publishing positive replication
Crneg
sigma.r = 
effort <- seq(1, 100, length.out = 1e4)
newHyp <- function(eta = 0.2, effort){
  h <- 1- eta*log10(effort)
  return(h)
}

pnh <- newHyp(effort = effort)


plot(effort, pnh, type = "l", ylab = "probability of new hypothesis")

plot(pnh, effort, type = "l", xlab = "probability of new hypothesis")
