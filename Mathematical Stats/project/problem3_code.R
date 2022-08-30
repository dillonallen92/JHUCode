# Problem 3, Glaucoma 

disease <- read.csv(file="glaucoma.csv", sep = ",", header = TRUE)
attach(disease)

result = c()

for(i in 1:62){
  x <- disease[,i][Class == "normal"]
  y <- disease[,i][Class == "glaucoma" ]
  result[i] <- t.test(x,y,alternative="two.sided")$p.value
}

pvalTable <- cbind(colnames(disease)[1:62], result)

names <- colnames(disease)
eyebox <- function(n){
  m <- n + 14
  par(mfrow=c(3,5))
  for(i in n:m){
    boxplot(disease[,i]~Class, main=names[i])
  }
}

eyebox(1)
eyebox(16)
eyebox(46)

# make the fit 
library(rpart)
fit <- rpart(Class~., data=disease, method="class")
fit

prediction <- predict(fit)
head(prediction)
gCorr  <- sum(prediction[,1] <= 0.5 & Class == 'glaucoma')
gIncor <- sum(prediction[,1] > 0.5  & Class == 'glaucoma')
nCorr  <- sum(prediction[,2] <= 0.5 & Class == 'normal')
nIncor <- sum(prediction[,2] > 0.5  & Class == 'normal')


