# Probabilities
pa = 0.5
pb = 0.3
pc = 0.2

# p that a won 2, b won 5, c won 3
na = 2
nb = 5
nc = 3
n = 10

pA <- choose(n,na)*((pa)^na)*(1-pa)^(n-na)
pB <- choose(n,nb)*((pb)^nb)*((1-pb)^(n-nb))
pC <- choose(n,nc)*((pc)^nc)*((1-pc)^(n-nc))
p <- pA + pB + pC
