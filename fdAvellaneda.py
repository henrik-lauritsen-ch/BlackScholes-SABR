
# import sys
# sys.path.append('../GitHub/BlackScholes-SABR')
# from BlackScholes import GarmanKohlhagen

import BlackScholes as bs
import numpy as np
import fdTools as fdT

 
class OptionPayOff:
    
    def __init__(self, optionType, exoticType) -> None:
        pass

    def PayOff(self):
        pass

class FiniteDifference2D:
    
    def __init__(self) -> None:
        pass


class FiniteDifferenceBS(FiniteDifference2D):
    
    def __init__(self) -> None:
        super().__init__()



option = bs.GarmanKohlhagen(100.0, 100.0, 1.0, 0.05, 0.100, 0.1)
ov = bs.Vanilla(100.0, 100.0, 1.0, 0.05, 0.10, 0.1)

######################################################################
# Caclulate dx, dt, maxS, minS
######################################################################
# t_N: number of steps in time, x_M: number of steps in space
x_interval = 360
x_M = x_interval + 1
t_N = 300 #later we choose time steps as function of space steps

# Todo: How to place Strike in grid  x_i < strike < x_i+1
# Todo: Somethig about move micro time step and smooth payoff in kink!? ... Tavella said something about this
# Todo: Proper calc max and min
max_spot = 160
min_spot = 40
delta_x = (max_spot - min_spot)/x_interval
delta_t = ov._expiryTerm/t_N


######################################################################
# Calculate Payoff Vector
######################################################################
option_type = bs.OptionType.Put
pay_off = np.empty(x_M, dtype=float)
x_vec = np.empty(x_M, dtype=float)

for i, payoff in enumerate(pay_off):
    spot_i = max_spot - i*delta_x
    x_vec[i] = spot_i
    # pay_off[i] = max(spot_i - ov._strike, 0)    
    pay_off[i] = max(ov._strike - spot_i, 0)    

print(x_vec, 'x_vec')
print(pay_off, 'pay_off t_0')


################################################################# 
#  Define functions a_h(x), b_h(x) and c_h(x)
################################################################# 
def a_h(a, b, i, x_vec):   
    x = x_vec[i]
    return -(a*x + b*x*x)
    
def c_h(a, b, i, x_vec):   
    x = x_vec[i]
    return (a*x - b*x*x)
    
def b_h(r_domestic, b, delta_t, i, x_vec):   
    x = x_vec[i]
    return (1 + r_domestic*delta_t + 2*b*x*x)


a = (ov._depositDomestic - ov._depositForeign)*delta_t/(2*delta_x)
b = 0.5*(ov._volatility*ov._volatility)*delta_t/(delta_x*delta_x)
r = np.array(pay_off)
print(r, 'r_0')

a_vec = np.empty(x_M - 1)
b_vec = np.empty(x_M)
c_vec = np.empty(x_M - 1)

# Build a, b, c
for i in range(0, x_M - 2):
    a_vec[i] = a_h(a, b, i+1, x_vec)

a_vec[x_M - 2] = -2*a*min_spot

for i in range(1, x_M - 1):
    b_vec[i] = b_h(ov._depositDomestic, b, delta_t, i, x_vec)
    c_vec[i] = c_h(a, b, i, x_vec)

b_vec[0] = 1 + ov._depositDomestic*delta_t - 2*a*max_spot
b_vec[x_M - 1] = 1 + ov._depositDomestic*delta_t + 2*a*min_spot
c_vec[0] = 2*a*max_spot
print(a_vec, 'a_vec')
print(b_vec, 'b_vec')
print(c_vec, 'c_vec')

########################################################
# Time Step
########################################################
for n in range(1, t_N + 1):   
    a = list(a_vec)
    b = list(b_vec)
    c = np.array(c_vec)    
    fdT.Thomas(a, b, c, r, x_M) #r: is passed by reference automatically so we can just run the function.
    
    

# a = np.identity(11)
# d = a@pay_off
# print(r)

import pandas as pd
file_location = '/Users/henriklauritsen/Documents/GitHub/BlackScholes-SABR/'
file_df = 'fd_test1.xlsx'


# df = pd.DataFrame(index = x_vec)

df = pd.DataFrame(x_vec, columns=['x_vec'])
df['pay_off'] = pay_off
df['r'] = r
df.to_excel(file_location + file_df)

print(ov.GetOptionValue(bs.OptionType.Put), 'BS call price')
print(delta_t/(delta_x*delta_x), 'dt/(dx*dx)')
print((max_spot - ov._strike)/delta_x, 'index')
print(df.loc[(max_spot - ov._strike)/delta_x, 'r'], 'Strike = 100')

# cubic spline prisen + läg grit omkring strike + lille timestep --> smooth price curve
# apply cvzone for text in count cars
# läs tavella kapitel 3