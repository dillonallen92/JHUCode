# golf game, find sigma^2 with a confidence interval for it
scores <- c(80, 88, 89, 90, 83, 84, 79, 82, 87, 76, 81, 86, 77, 78, 85)

chisq_lb <- 26.119
chisq_ub <- 5.629
var_scores <- var(scores)
n = length(scores)

lb = (n-1)*var_scores/chisq_lb
ub = (n-1)*var_scores/chisq_ub
