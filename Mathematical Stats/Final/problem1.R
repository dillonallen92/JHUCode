# problem 1

ants <- rep(0:5, times=c(10,29,25,27,6,3))
table(ants)

# MLE for poisson is mean
lambda = mean(ants)
n = 100

# estimated values
est_arr = c()
for (i in 0:5){
  est_arr[i+1] = n*exp(-lambda)*lambda^(i)/factorial(i)
}

sumprobs = 0
for( k in 0:5){
  sumprobs = sumprobs + exp(-lambda)*lambda^k / factorial(k)
}

est_arr[7] = n*(1-sumprobs)

obs = c(10,29,25,27,6,3)
est_arr_redone = c(est_arr[1], est_arr[2], est_arr[3], est_arr[4], est_arr[5], est_arr[6] + est_arr[7])

d = 0
for(i in 1:6){
  d = d + (obs[i] - est_arr_redone[i])^2/obs[i]
}

