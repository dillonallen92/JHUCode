# Quality Control Stuff
A <- c(31, 45, 52, 30)
B <- c(62, 35, 30, 40)
C <- c(20, 20, 25, 23)

t.test(x=A, y=B, paired = TRUE)
t.test(x=A, y=C, paired = TRUE)
t.test(x=B, y=C, paired = TRUE)
