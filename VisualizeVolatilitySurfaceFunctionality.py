

# 2. SABR Calibration + find solver!!
# 3. Print smile med Cubic-spline + SABR calibration under main()
# 4. Pricing Vanilla Options Under SABR model = Combine BS + 2
# 5. Pricing Exotic Options under SABR ... study!!


"""
 FX Vanilla option tools: Show impact of BS with cubic spline and under SABR interpolation

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
expiryTerm = 720/365
r = 0.002978
q = 0.007450
spot = 76.1340
vol_smile = np.array([0.117885, 0.1191, 0.1300, 0.1501, 0.174995])
beta = 0.85
# vol_smile = np.array([0.14852, 0.10042, 0.0973, 0.11582, 0.23732])
strike = 75.96442


                      
# Strike from Delta points
sd = sfd.StrikeFromDelta(spot, r, q, expiryTerm)
strike_vec = sd.GetStrikeVector(vol_smile)
# logmon_strike_vec = sd.GetLogMoneynessStrikeVector(vol_smile)

# SABR first calibration
sabr = vs.SABRVolSurface(spot, r, q, expiryTerm, vol_smile)
# sabr.SabrCalibration()
print(sabr.GetVolatility(strike), 'test')

vovol = sabr._vovol
alpha =  sabr._alpha
corr = sabr._corr
print('alpha0: ', sabr._alpha0, ' corr0: ', sabr._corr0, ' vovol0: ', sabr._vovol0, ' beta ', sabr._beta)

print('alpha: ', sabr._alpha, ' corr: ', sabr._corr, ' vovol: ', sabr._vovol, ' beta ', sabr._beta)
# print(sabr._alpha, sabr._corr, sabr._vovol, sabr._beta)



# BlackScholes Option object
atm_strike = sd.GetATMStrike(vol_smile[2])
bs_option = bs.Vanilla(spot, atm_strike, expiryTerm, r, q, vol_smile[2])



option_price_vec = np.zeros(5)
for i in range(5):
    bs_option._strike = strike_vec[i] #math.exp(logmon_strike_vec[i])*atm_strike
    option_price_vec[i] = bs_option.GetBaseOptionValue(bs.OptionType.Call, vol_smile[i])    
    

plot_strikes = np.linspace(strike_vec[0]*0.975, strike_vec[4]*1.025, 60)
dt = pd.DataFrame(plot_strikes, columns=['Strikes'])
cs = u.CubicSplineInterpolation()
pl = u.PiecewiseLinearInterpolation()

# dt['Actual_Strike'] = dt.apply(lambda row: math.exp(row['Strikes'])*atm_strike, axis=1) 
dt['CS_Vol_Smile'] = dt.apply(lambda row: cs.GetInterpolatedValue(row['Strikes'], strike_vec, vol_smile), axis=1)
dt['PL_Vol_Smile'] = dt.apply(lambda row: pl.GetInterpolatedValue(row['Strikes'], strike_vec, vol_smile), axis=1)
# dt['SABR_Vol_FirstG'] = dt.apply(lambda row: sabr.SabrImpliedVol(row['Strikes'], sabr._corr0, sabr._vovol0, sabr._alpha0, 0.85), axis=1)
dt['SABR_Vol_Smile'] = dt.apply(lambda row: sabr.SabrImpliedVol(row['Strikes'], corr, vovol, alpha, 0.85), axis=1)

# corr *= 0.8
# dt['SC_SABR_Vol_Smile'] = dt.apply(lambda row: sabr.SabrImpliedVol(row['Strikes'], corr, vovol, alpha, 0.85), axis=1)

for index, i in dt.iterrows():        
    bs_option._strike = math.exp(dt.loc[index, 'Strikes'])*atm_strike
    dt.loc[index, 'OptionPrice'] = bs_option.GetBaseOptionValue(bs.OptionType.Call, dt.loc[index, 'CS_Vol_Smile'])
    dt.loc[index, 'PL_OptionPrice'] = bs_option.GetBaseOptionValue(bs.OptionType.Call, dt.loc[index, 'PL_Vol_Smile'])
    dt.loc[index, 'InnerValue'] = max(bs_option._spot - bs_option._strike, 0)


# print(dt)

fig = plt.figure(figsize=(11, 9.5))
fig.subplots_adjust(hspace=0.275)

ax =fig.add_subplot(111)
# ax =fig.add_subplot(211)
ax.plot(dt['Strikes'], dt['CS_Vol_Smile'], label = 'Cubic Spline Smile')
ax.plot(dt['Strikes'], dt['PL_Vol_Smile'], label = 'Piecewise Linear Smile')
ax.plot(dt['Strikes'], dt['SABR_Vol_Smile'], label = 'SABR calib')
# ax.plot(dt['Strikes'], dt['SABR_Vol_FirstG'], label = 'SABR first guess')
# ax.plot(dt['Strikes'], dt['SC_SABR_Vol_Smile'], label = 'SABR manual corr')
ax.scatter(strike_vec, vol_smile, facecolor = 'orange', linewidth=2.1, label = 'Quoted Points')
ax.set_title('EURNOK, Vol-Smile')
ax.set_xlabel('Strike')
ax.set_ylabel('Volatility')
ax.legend()
ax.grid(True)


# ax = fig.add_subplot(212)
# ax.plot(dt['Strikes'], dt['OptionPrice'], label = 'Cubic Spline Prices')
# ax.plot(dt['Strikes'], dt['PL_OptionPrice'], label = 'Piecewise Linear Prices')
# # ax.plot(dt['Strikes'], dt['InnerValue'], label = 'Option: Inner Value')
# ax.scatter(strike_vec, option_price_vec, color='orange', linewidth=2.1, label='Quoted Points')
# ax.set_title('EURNOK Call Price along the volatility smile')
# ax.set_xlabel('Log-Moneyness')
# ax.set_ylabel('Option Price (NOK/Option)')
# ax.legend()
# ax.grid(True)

plt.show()

# print(dt)


strikes = np.array([2, 4, 6, 8, 10, 12, 14])
dtf = pd.DataFrame(strikes, columns=['Strikes'])
sabr2 = vs.SABRVolSurface(8.01, 0, 0, 15.0, np.array([0.175, 0.08, 0.05, 0.06, 0.07]), 0.4)

dtf['SABR2_Vol_Smile'] = dt.apply(lambda row: sabr2.SabrImpliedVol(row['Strikes'], -0.33, 0.25, 0.0825, 0.4), axis=1)
print(dtf)
# BELOW FOR TESTING PURPOSES!!!


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