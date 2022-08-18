fish <- read.csv("fish.csv",sep = ",",header = TRUE)

# pair test
t.test(fish$Guarding, fish$Fanning, paired = TRUE)

# sign test
guarding <- fish$Guarding
fanning <- fish$Fanning
diff <- fish$Guarding - fish$Fanning
