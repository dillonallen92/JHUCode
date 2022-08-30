# two sample mean hypothesis test

C <- c(4.62, 5.95, 4.81, 5.26, 5.41, 5.64, 4.84, 3.96, 6.38, 5.24, 5.07, 3.63)
CM <- c(3.86, 10.95, 5.85, 4.61, 7.90, 6.31, 7.07, 5.72, 4.60, 6.66, 5.73, 6.97, 4.26, 3.62, 4.71)

# alpha = 0.05 level
# H0 = muCM = muC
# HA = muCM > muC

n = length(CM)
m = length(C)
x = mean(CM)
y = mean(C)
varx <- var(CM)
vary <- var(C)
sp = sqrt(((n-1)*varx + (m-1)*vary)/(n+m-2))
df = n + m - 2
t = (x-y)/(sp*(sqrt((1/n) + (1/m))))
