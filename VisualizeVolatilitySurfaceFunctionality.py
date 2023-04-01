
"""
 FX Vanilla option tools: Show impact of BS with cubic spline and under SABR interpolation

HISTORIC TEST DATA:
CHF/DKK	Spot	5.7444					
                                                               Dom ccy DKK	For ccy CHF
7   vol_smile = np.array([0.120, 0.115, 0.113, 0.114, 0.118])  0.006056202	0.000547894
30  vol_smile = np.array([0.120, 0.116, 0.115, 0.118, 0.124])  0.007074041	0.001127818
181 vol_smile = np.array([0.112, 0.112, 0.115, 0.125, 0.137])  0.011553891	0.004381738
365 vol_smile = np.array([0.109, 0.109, 0.115, 0.128, 0.147])  0.015856035	0.00669313

ZAR/JPY	Spot	11.7336						
                                                                Dom ccy JPY	 For ccy ZAR
7   vol_smile = np.array([0.188, 0.172, 0.155, 0.148, 0.148])  0.000901339	0.056517435
30  vol_smile = np.array([0.189, 0.171, 0.153, 0.146, 0.146])  0.001671885	0.061701683
181 vol_smile = np.array([0.233, 0.200, 0.172, 0.161, 0.158])  0.003722718	0.061610374
365 vol_smile = np.array([0.250, 0.209, 0.179, 0.171, 0.169])  0.00503213	0.064996107
 
 USDJPY	Spot	82.82						
                                                                Dom ccy JPY	 For ccy USD
7   vol_smile = np.array([0.110, 0.103, 0.101, 0.106, 0.114])  0.000901339	0.003041578
31  vol_smile = np.array([0.114, 0.106,	0.101, 0.103, 0.109])  0.001671885	0.004316413
179 vol_smile = np.array([0.139, 0.124, 0.116, 0.115, 0.121])  0.003722718	0.007776803
365 vol_smile = np.array([0.154, 0.136, 0.126, 0.123, 0.131])  0.00503213	0.010589565


"""

import pandas as pd
import numpy as np
import BlackScholes as bs
import matplotlib.pyplot as plt
import VolatilitySurface as vs
import Utility as u
import math
import StrikeFromDelta as sfd

# Volatility Smile
# spot = 76.1340
# strike = 75.96442
# expiryTerm = 365/365
# r = 0.002978
# q = 0.007450
beta = 0.85
# vol_smile = np.array([0.117885, 0.1191, 0.1300, 0.1501, 0.174995])
# vol_smile = np.array([0.14852, 0.10042, 0.0973, 0.11582, 0.23732])
# vol_smile = np.array([0.09852, 0.09542, 0.0973,	0.10582, 0.11732])

# USDJPY
spot = 11.7336
strike = 11.7336
expiryTerm = 7/365.0
r = 0.000901339
q = 0.056517435

#7, 30, 181, 365
vol_smile = np.array([0.188, 0.172, 0.155, 0.148, 0.148])  #0.000901339	0.056517435
# vol_smile = np.array([0.189, 0.171, 0.153, 0.146, 0.146])  #0.001671885	0.061701683
# vol_smile = np.array([0.233, 0.200, 0.172, 0.161, 0.158])  #0.003722718	0.061610374
# vol_smile = np.array([0.250, 0.209, 0.179, 0.171, 0.169])  #0.00503213	0.064996107

                      
# Strike from Delta Object
sd = sfd.StrikeFromDelta(spot, r, q, expiryTerm)
strike_vec = sd.GetStrikeVector(vol_smile)
atm_strike = strike_vec[2]

# SABR-Surface Object
sabr = vs.SABRVolSurface(spot, r, q, expiryTerm, vol_smile)
sabr.SabrCalibration()
print('alpha0: ', sabr._alpha0, ' corr0: ', sabr._corr0, ' vovol0: ', sabr._vovol0, ' beta ', sabr._beta)
print('alpha: ', sabr._alpha, ' corr: ', sabr._corr, ' vovol: ', sabr._vovol, ' beta ', sabr._beta)


# BlackScholes object
bs_option = bs.Vanilla(spot, atm_strike, expiryTerm, r, q, vol_smile[2])
option_price_vec = np.zeros(5)
for i in range(5):
    bs_option._strike = strike_vec[i] #math.exp(logmon_strike_vec[i])*atm_strike
    option_price_vec[i] = bs_option.GetBaseOptionValue(bs.OptionType.Call, vol_smile[i])    
    

# Granular x-axis for plot
plot_strikes = np.linspace(strike_vec[0]*0.975, strike_vec[4]*1.025, 60)
plot_moneyness_vec = plot_strikes/atm_strike*100


dt = pd.DataFrame(plot_strikes, columns=['Strikes'])
dt['Moneyness'] = plot_moneyness_vec
cs = u.CubicSplineInterpolation()
pl = u.PiecewiseLinearInterpolation()

# dt['Actual_Strike'] = dt.apply(lambda row: math.exp(row['Strikes'])*atm_strike, axis=1) 
dt['CS_Vol_Smile'] = dt.apply(lambda row: cs.GetInterpolatedValue(row['Strikes'], strike_vec, vol_smile), axis=1)
dt['PL_Vol_Smile'] = dt.apply(lambda row: pl.GetInterpolatedValue(row['Strikes'], strike_vec, vol_smile), axis=1)
dt['SABR_Vol_Smile'] = dt.apply(lambda row: sabr.SabrImpliedVol(row['Strikes'], sabr._alpha, sabr._corr, sabr._vovol, 0.85), axis=1)


for index, i in dt.iterrows():        
    bs_option._strike = dt.loc[index, 'Strikes']
    dt.loc[index, 'OptionPrice'] = bs_option.GetBaseOptionValue(bs.OptionType.Call, dt.loc[index, 'CS_Vol_Smile'])
    dt.loc[index, 'PL_OptionPrice'] = bs_option.GetBaseOptionValue(bs.OptionType.Call, dt.loc[index, 'PL_Vol_Smile'])
    dt.loc[index, 'SABR_OptionPrice'] = bs_option.GetBaseOptionValue(bs.OptionType.Call, dt.loc[index, 'SABR_Vol_Smile'])
    dt.loc[index, 'InnerValue'] = max(bs_option._spot - bs_option._strike, 0)

print(dt)

fig = plt.figure(figsize=(11, 9.5))
fig.subplots_adjust(hspace=0.275)
strike_vec = strike_vec/atm_strike*100


ax =fig.add_subplot(111)
# ax =fig.add_subplot(211)
ax.plot(dt['Moneyness'], dt['CS_Vol_Smile'], label = 'Cubic spline Smile')
ax.plot(dt['Moneyness'], dt['PL_Vol_Smile'], label = 'Piecewise linear Smile')
ax.plot(dt['Moneyness'], dt['SABR_Vol_Smile'], label = 'SABR Smile')
ax.scatter(strike_vec, vol_smile, facecolor = 'orange', linewidth=2.1, label = 'Quoted Points')
ax.set_title('ZARJPY, 1W Volatility Smile')
ax.set_xlabel('Moneyness')
ax.set_ylabel('Volatility')
ax.legend()
ax.grid(True)


# ax = fig.add_subplot(212)
# ax.plot(dt['Log-Moneyness'], dt['OptionPrice'], label = 'Cubic Spline Prices')
# ax.plot(dt['Log-Moneyness'], dt['PL_OptionPrice'], label = 'Piecewise Linear Prices')
# ax.plot(dt['Log-Moneyness'], dt['SABR_OptionPrice'], label = 'SABR Price')
# ax.scatter(strike_vec, option_price_vec, color='orange', linewidth=2.1, label='Quoted Points')
# ax.set_title('EURNOK Call Price along the volatility smile')
# ax.set_xlabel('Log-Moneyness')
# ax.set_ylabel('Option Price (NOK/Option)')
# ax.legend()
# ax.grid(True)

plt.show()





# BELOW FOR TESTING PURPOSES!!!

# print(dt)
# strikes = np.array([2, 4, 6, 8, 10, 12, 14])
# dtf = pd.DataFrame(strikes, columns=['Strikes'])
# sabr2 = vs.SABRVolSurface(8.01, 0, 0, 15.0, np.array([0.175, 0.08, 0.05, 0.06, 0.07]), 0.4)

# dtf['SABR2_Vol_Smile'] = dt.apply(lambda row: sabr2.SabrImpliedVol(row['Strikes'], -0.33, 0.25, 0.0825, 0.4), axis=1)
# print(dtf)



#print(dt)
# print(plot_vol)
# ss = np.arange(0.1, 1.0, 0.05)
# print(ss)
# plot_vol = np.array(range(30))
# option = bs.Vanilla(100, 95, 0.5, 0.01, 0.02, 0.30)
# print(option.GetOptionValue(bs.OptionType.Call))

# vol_x = np.linspace(0.025, 1.425, 60)
# df = pd.DataFrame(index=vol_x, columns=range(1))

# df.columns = ['Volatility']
# vol = np.linspace(0.025, 1.425, 60)
# df['Volatility'] = vol

# for index, row in df.iterrows():
#     option._volatility = index #df.loc[index, 'Volatility']
#     df.loc[index, 'OptionValue'] = option.GetOptionValue(bs.OptionType.Call)
#     df.loc[index, 'Vega'] = option.GetVega()
#     df.loc[index, 'Volga'] = option.GetVolga()

#   for index, row in data.iterrows():
#         data.loc[index,('hold_index_' + strat_nr)] = index_yesterday


# df_ir_hs = pd.DataFrame(index=range(num_ir_hs), columns=range(6))
# df_ir_hs.columns = ['Hedge_Set_Id', 'IR_Currency', 'Market_Value', 'Collateral_Cover_Value','AdjustedNotional', 'Add_On_HS']






# 5. Pricing Exotic Options under SABR ... study!!
