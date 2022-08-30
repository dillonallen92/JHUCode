# Finding the MLE for some distribution

n = 10
sum1 = 0
sum2 = 0
sum3 = 0
sum4 = 0
for (i in 1:n){
  vals1 <- runif(10, min=0, max=1) * (2/3)
  vals2 <- runif(100, min=0, max=1) * (2/3)
  vals3 <- runif(1000, min = 0, max = 1) * (2/3)
  vals4 <- runif(10000, min = 0, max = 1) * (2/3)
  
  min1 <- min(vals1)
  min2 <- min(vals2)
  min3 <- min(vals3)
  min4 <- min(vals4)
  
  sum1 = sum1 + min1
  sum2 = sum2 + min2
  sum3 = sum3 + min3
  sum4 = sum4 + min4
}

avg1 = sum1/n
avg2 = sum2/n
avg3 = sum3/n
avg4 = sum4/n