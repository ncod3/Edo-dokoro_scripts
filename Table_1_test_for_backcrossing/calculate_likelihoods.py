#! /usr/bin/env python3

import math


AA = 434253
Aa = 460070
aa = 20760
p  = 0.02


C1 = math.log(math.factorial(AA + Aa + aa))
C2 = math.log(math.factorial(AA)) + math.log(math.factorial(Aa)) + math.log(math.factorial(aa))
C = C1 - C2

BC1F1_lnl = C + AA*math.log((1- p)/2) + Aa*math.log((1- p)/2) + aa*math.log(p)
BC2F1_lnl = C + AA*math.log(3*(1- p)/4) + Aa*math.log((1 - p)/4) + aa*math.log(p)

print(BC1F1_lnl)
print(BC2F1_lnl)
