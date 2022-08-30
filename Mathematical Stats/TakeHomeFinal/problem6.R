x <- c(1.00, 6.00,  9.00, 12.00, 14.00, 15.0, 17.00, 18.00, 19.00, 20.00)
y <- c(3.06, 7.23, 19.34, 29.29, 30.76, 36.5, 39.24, 40.54, 52.19, 46.98)

# (a) Fit to a simple linear model
# Use theorem 11.3.1 for b0 (intercept) and b1 (slope)

n = length(x)
sumxy = 0
sumx = sum(x)
sumy = sum(y)
sumx2 = 0
sumy2 = 0
sumxx2 = 0
for(i in 1:n){
  sumxy = sumxy + x[i]*y[i]
  sumx2 = sumx2 + x[i]*x[i]
  sumy2 = sumy2 + y[i]*y[i]
  sumxx2 = sumxx2 + (x[i] - mean(x))^2
}

b1 = (n*sumxy - (sumx)*(sumy))/(n*sumx2 - (sumx^2))
b0 = mean(y) - b1*mean(x)

# (b) Conduct t-test at alpha = 0.05 for slope being equal to 0
# H0: b1 = 0
# Ha: b1 != 0

b1prime = 0
s2 = (1/(n-2))*(sumy2 - b0*sumy - b1*sumxy)
s = sqrt(s2)
t = (b1 - b1prime)/(s/sqrt(sumxx2))
# t_(alpha/2,8) = +/- 2.3060
t_alpha = 2.3060
# since t > t_(alpha/2,8), we can reject the null hypothesis that the slope was 0

# (c) develop a 95% confidence interval for the intercept 
width = t_alpha*((s * sqrt(sumx2))/(sqrt(n)*sqrt(sumxx2)))

lb = b0 - width
ub = b0 + width

# (d) develop a confidence interval for the regression line when x = 5
x_val = 5
y_hat = b0 + b1*x_val

width_rl = t_alpha * s * sqrt((1/n) + (x_val - mean(x))/sumxx2)

lb_rl = y_hat - width_rl
ub_rl = y_hat + width_rl

# (e) 

width_newval = t_alpha * s * sqrt(1 + (1/n) + (x_val - mean(x))/sumxx2)

lb_newval = y_hat - width_newval
ub_newval = y_hat + width_newval