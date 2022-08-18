fish <- read.csv("fish.csv",sep = ",",header = TRUE)

# pair test
t.test(fish$Guarding, fish$Fanning, paired = TRUE)

# sign test
guarding <- fish$Guarding
fanning <- fish$Fanning
diff <- fish$Guarding - fish$Fanning
# choosing an alpha of 0.05
z_alpha = 1.96
n = length(diff)
k = sum(diff > 0)
z = (k - n/2)/(sqrt(n/4))

# split into groups now
fish_split <- split(fish, fish$pH)
fish_ph55 <- fish_split$`5.5`
fish_ph6 <- fish_split$`6`
fish_ph65 <- fish_split$`6.5`
fish_ph7 <- fish_split$`7`

# tests for pH 5.5 
t.test(fish_ph55$Guarding, fish_ph55$Fanning, paired = TRUE) # t = 9.8986, alt hyp
binom.test(sum((fish_ph55$Guarding - fish_ph55$Fanning) > 0), length(fish_ph55$Guarding)) #p-val = 4.005e-5

# tests for ph6
t.test(fish_ph6$Guarding, fish_ph6$Fanning, paired = TRUE) # t = 4.1506, alt hyp
binom.test(sum((fish_ph6$Guarding - fish_ph6$Fanning) > 0), length(fish_ph6$Guarding)) #p-val = 0.04139

# test for ph6.5
t.test(fish_ph65$Guarding, fish_ph65$Fanning, paired = TRUE) # t = 5.3816, alt hyp
binom.test(sum((fish_ph65$Guarding - fish_ph65$Fanning) > 0), length(fish_ph65$Guarding)) #p-val = .0004025

# test for ph7
t.test(fish_ph7$Guarding, fish_ph7$Fanning, paired = TRUE) # t = 5.3816, alt hyp
binom.test(sum((fish_ph7$Guarding - fish_ph7$Fanning) > 0), length(fish_ph7$Guarding)) #p-val = .8238

# Split by genders
fish_gender_split <- split(fish, fish$Gender)
fish_male <- fish_gender_split$M
fish_female <- fish_gender_split$F

# male test
t.test(fish_male$Guarding, fish_male$Fanning, paired = TRUE) # p-val = 2.2e-16
binom.test(sum((fish_male$Guarding - fish_male$Fanning)>0), length(fish_male$Guarding)) #p-val = 1.8e-12

# female test
t.test(fish_female$Guarding, fish_female$Fanning, paired=TRUE) #p-val = 0.1455
binom.test(sum((fish_female$Guarding - fish_female$Fanning)>0), length(fish_female$Guarding)) #p-val = 0.4296

# now it is time to start plotting things
library(ggplot2)













