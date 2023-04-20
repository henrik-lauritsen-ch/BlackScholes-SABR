
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
import unittest
import pandas as pd
import numpy as np
import BlackScholes as bs
import matplotlib.pyplot as plt
import VolatilitySurface as vs
import SABRWing as sw
import Utility as u
# import math
import StrikeFromDelta as sfd

beta = 0.85
spot = 5.7444
strike = 5.7444
expiryTerm = 30/365.0
r = 0.007074041
q = 0.001127818

# 7, 30, 181, 365                                              Dom ccy DKK	For ccy CHF
# vol_smile = np.array([0.120, 0.115, 0.113, 0.114, 0.118])  # 0.006056202	0.000547894
# vol_smile = np.array([0.120, 0.116, 0.115, 0.118, 0.124])  # 0.007074041	0.001127818
# vol_smile = np.array([0.112, 0.112, 0.115, 0.125, 0.137])  # 0.011553891	0.004381738
# vol_smile = np.array([0.109, 0.109, 0.115, 0.128, 0.147])  # 0.015856035	0.00669313
vol_smile = np.array([0.250, 0.209, 0.179, 0.171, 0.169])
                      
# Strike from Delta Object
sd = sfd.StrikeFromDelta(spot, r, q, expiryTerm)
strike_vec = sd.GetStrikeVector(vol_smile)
atm_strike = strike_vec[2]

# SABR-Surface Object
sabr = vs.SABRVolSurface(spot, r, q, expiryTerm, vol_smile)
sabr_wing = sw.SABRWingSurface(spot, r, q, expiryTerm, vol_smile)

# BlackScholes Object
bs_option = bs.Vanilla(spot, atm_strike, expiryTerm, r, q, vol_smile[2])


option_price_vec = np.zeros(5)
for i in range(5):
    bs_option._strike = strike_vec[i] #math.exp(logmon_strike_vec[i])*atm_strike
    option_price_vec[i] = bs_option.GetBaseOptionValue(bs.OptionType.Call, vol_smile[i])    


# Granular x-axis for plot
plot_strikes = np.linspace(strike_vec[0]*0.95, strike_vec[4]*1.05, 120)
delta_K = plot_strikes[1] - plot_strikes[0]
plot_moneyness_vec = plot_strikes/atm_strike*100

# Create DataFrame
dt = pd.DataFrame(plot_strikes, columns=['Strikes'])
dt['Moneyness'] = plot_moneyness_vec
cs = u.CubicSplineInterpolation()
pl = u.PiecewiseLinearInterpolation()


dt['CS_Vol_Smile'] = dt.apply(lambda row: cs.GetInterpolatedValue(row['Strikes'], strike_vec, vol_smile), axis=1)
dt['PL_Vol_Smile'] = dt.apply(lambda row: pl.GetInterpolatedValue(row['Strikes'], strike_vec, vol_smile), axis=1)
dt['SABR_Vol_Smile'] = dt.apply(lambda row: sabr.GetVolatility(row['Strikes']), axis=1)
dt['SABR_Wing_Vol'] = dt.apply(lambda row: sabr_wing.GetVolatility(row['Strikes']), axis=1)


for index, i in dt.iterrows():        
    bs_option._strike = dt.loc[index, 'Strikes']
    dt.loc[index, 'CS_Price'] = bs_option.GetBaseOptionValue(bs.OptionType.Call, dt.loc[index, 'CS_Vol_Smile'])
    dt.loc[index, 'PL_Price'] = bs_option.GetBaseOptionValue(bs.OptionType.Call, dt.loc[index, 'PL_Vol_Smile'])
    dt.loc[index, 'SABR_Price'] = bs_option.GetBaseOptionValue(bs.OptionType.Call, dt.loc[index, 'SABR_Vol_Smile'])
    dt.loc[index, 'SABR_WingPrice'] = bs_option.GetBaseOptionValue(bs.OptionType.Call, dt.loc[index, 'SABR_Wing_Vol'])
    # dt.loc[index, 'InnerValue'] = max(bs_option._spot - bs_option._strike, 0)


for index, i in dt.iterrows():
    
    if ((index >= 1) and (index <= 119 - 1)):
        dt.loc[index, 'CS_pdf'] = (dt.loc[index + 1, 'CS_Price']+ dt.loc[index - 1, 'CS_Price'] - 2*dt.loc[index, 'CS_Price'])/(delta_K*delta_K)
        dt.loc[index, 'PL_pdf'] = (dt.loc[index + 1, 'PL_Price']+ dt.loc[index - 1, 'PL_Price'] - 2*dt.loc[index, 'PL_Price'])/(delta_K*delta_K)
        dt.loc[index, 'SABR_pdf'] = (dt.loc[index + 1, 'SABR_Price']+ dt.loc[index - 1, 'SABR_Price'] - 2*dt.loc[index, 'SABR_Price'])/(delta_K*delta_K)
        dt.loc[index, 'SABRWing_pdf'] = (dt.loc[index + 1, 'SABR_WingPrice']+ dt.loc[index - 1, 'SABR_WingPrice'] - 2*dt.loc[index, 'SABR_WingPrice'])/(delta_K*delta_K)

# print(dt)


fig = plt.figure(figsize=(11, 9.5))
fig.subplots_adjust(hspace=0.275)
strike_vec = strike_vec/atm_strike*100

# ax =fig.add_subplot(111)
# ax.plot(dt['Moneyness'], dt['CS_Vol_Smile'], label = 'Cubic spline Smile')
# # ax.plot(dt['Moneyness'], dt['PL_Vol_Smile'], label = 'Piecewise linear Smile')
# ax.plot(dt['Moneyness'], dt['SABR_Vol_Smile'], label = 'SABR Smile')
# ax.plot(dt['Moneyness'], dt['SABR_Wing_Vol'], label = 'SABR Wing')
# ax.scatter(strike_vec, vol_smile, facecolor = 'orange', linewidth=2.1, label = 'Quoted Points')
# ax.set_title('CHF/DKK, 1-Month Volatility Smile')
# ax.set_xlabel('Moneyness')
# ax.set_ylabel('Volatility')
# ax.legend()
# ax.grid(True)


ax = fig.add_subplot(111)
ax.plot(dt['Moneyness'], dt['CS_pdf'], label = 'Cubic Spline pdf')
# ax.plot(dt['Moneyness'], dt['PL_pdf'], label = 'Piecewise Linear pdf')
ax.plot(dt['Moneyness'], dt['SABR_pdf'], label = 'SABR pdf')
ax.plot(dt['Moneyness'], dt['SABRWing_pdf'], label = 'SABR Wing pdf')
# ax.scatter(strike_vec, option_price_vec, color='orange', linewidth=2.1, label='Quoted Points')
ax.set_title('CHF/DKK implied risk neutral distribution')
ax.set_xlabel('Moneyness')
ax.set_ylabel('pdf-function')
ax.legend()
ax.grid(True)

plt.show()


