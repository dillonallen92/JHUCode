# traffic circle debate
n = 182
accidents_before = 32 # acc_B
accidents_after = 20 # acc_A

# alpha = 0.05

# H0 pB = pA
# HA pB > pA

pB = accidents_before / n
pA = accidents_after / n
pE = (accidents_before + accidents_after)/(n + n)

# need to get a z > 1.96 for alpha = 0.05
z_num = pB - pA
z_denom = sqrt((pE*(1-pE))/n + (pE*(1-pE))/n)
z = z_num / z_denom