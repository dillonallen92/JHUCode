# viscosity problems

v_b <- c(20.65, 18.50, 20.08, 20.87, 19.46, 19.79, 20.01, 19.17, 20.26, 20.44, 19.38, 20.78)
v_c <- c(20.16, 18.92, 20.11, 20.32, 20.34, 19.85, 20.07, 20.14, 20.05, 20.09)

sum_b = sum(v_b)
sum_c = sum(v_c)
sum_b2 = 0
sum_c2 = 0
n = length(v_b)
m = length(v_c)
for(i in 1:n){
  sum_b2 = sum_b2 + v_b[i]^2
}

for(i in 1:m){
  sum_c2 = sum_c2 + v_c[i]^2
}

sb2 = ((n*sum_b2) - sum_b^2)/(n*(n-1))
sc2 = ((m*sum_c2) - sum_c^2)/(m*(m-1))

f = sc2/sb2
