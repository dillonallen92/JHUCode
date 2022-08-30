# Marksman shots

mux = 0.5
muy = 1
varx = 4
vary = 6
rho = 0.7
x = 1
d = 5 # diameter of circle

# equation of circle, x^2 + y^2 = r^2, find the upper and lower coordinates
cu = sqrt((d/2)^2 - 1)
cl = - sqrt((d/2)^2 - 1)

expy_givenx = muy + (rho*sqrt(vary))*(1 - mux)/sqrt(varx)
vary_givenx = (1-rho^2)*vary

# now find the z values
zl = (cl - expy_givenx)/sqrt(vary_givenx)
zu = (cu - expy_givenx)/sqrt(vary_givenx)

# P(cl < x < cu) = P(zu) - P(zl)
p = .3121 - .0170
