# Problem 2, Generic Problem set with Regression

generic <- read.csv("generic.csv", sep = ",", header = TRUE)

# getting tired of typing generic$(blah), substituting 
y <- generic$y
x1 <- generic$x1 
x2 <- generic$x2 
x3 <- generic$x3
x4 <- generic$x4
x5 <- generic$x5
x6 <- generic$x6
x7 <- generic$x7 

# exploratory data analysis 
par(mfrow=c(2,4), mai=c(0.3,0.1,0.2,0.1))
hist(generic$y, main = "Response Histogram") # right-skewed
hist(generic$x1, main = "x1 histogram") # about normal
hist(generic$x2, main = "x2 histogram") # about normal
hist(generic$x3, main = "x3 histogram") # 'uniform' but not
hist(generic$x4, main = "x4 histogram") # about normal
hist(generic$x5, main = "x5 histogram") # 'uniform'
hist(generic$x6, main = "x6 histogram") # almost bimodal?
hist(generic$x7, main = "x7 histogram") # straight up 0 or 1, maybe categorical

# Check some correlations
cor(x1, y) # very low, 0.03185
cor(x2, y) # negative, -0.4132
cor(x3, y) # positive, 0.41152
cor(x4, y) # small negative, -0.00868
cor(x5, y) # positive, big 0.5486
cor(x6, y) # low positive, 0.02811
cor(x7, y) # positive, 0.31876

# do some models and grab summaries from it
model1 = lm(y~x1)
summary(model1) # r^2 horrible
model2 = lm(y~x2)
summary(model2) # r^2 = .1623, ehh
model3 = lm(y~x3)
summary(model3) # r^2 = .1609
model4 = lm(y~x4)
summary(model4) # r2 = -0.01013, bad
model5 = lm(y~x5)
summary(model5) # r2 = 0.2938
model6 = lm(y~x6)
summary(model6) # r2 = -0.009405
model7 = lm(y~x7)
summary(model7) # r2 = 0.09244

# combined model
combinedModel <- lm(y~x1+x2+x3+x4+x5+x6+x7)
summary(combinedModel)

# the t-values that are most significant are x2 (-1.528), x3 (5.863),
# x5 (6.819), x7 (5.811)

res <- resid(combinedModel)

hist(res)

# plot the residuals vs each predictor variable
plot(x1, res)
plot(x2, res)
plot(x3, res)
plot(x4, res)
plot(x5, res)
plot(x6, res)
plot(x7, res)
# do some sort of transformation and check the model summary on that
y <- log(y)
testModel1 <- lm(y~ x2^2 + x3 + x5 + x7)
summary(testModel1)



