
# 1. Do piece wise - lin interpol and show how greeks move
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

# Volatility Smile
vol_smile = np.array([0.09852, 0.09542, 0.0973, 0.10582, 0.11732])

# Strike from Delta points
sfd = vs.StrikeFromDelta(10.3719, 0.00565, 0.01822, 71.0/365.0)
strike_vec = sfd.GetStrikeVector(vol_smile)

# BlackScholes Option object
bs_option = bs.Vanilla(10.3719, 10.3719, 71/365.0, 0.00565, 0.01822, 0.0973)

option_price_vec = np.zeros(5)
for i in range(5):
    bs_option._strike = strike_vec[i]
    option_price_vec[i] = bs_option.GetBaseOptionValue(bs.OptionType.Call, vol_smile[i])    
    



plot_strikes = np.linspace(9.75, 11.25, 30)
plot_vol = np.zeros(30)
dt = pd.DataFrame(plot_strikes, columns=['Strikes'])
dt['Vol_Smile'] = dt.apply(lambda row: u.CubicSplineInterpolation(row['Strikes'], strike_vec, vol_smile), axis=1)

for index, i in dt.iterrows():    
    bs_option._strike = dt.loc[index, 'Strikes']
    dt.loc[index, 'OptionPrice'] = bs_option.GetBaseOptionValue(bs.OptionType.Call, dt.loc[index, 'Vol_Smile'])
    dt.loc[index, 'InnerValue'] = max(bs_option._spot - bs_option._strike, 0)

fig = plt.figure(figsize=(11, 9.5))
fig.subplots_adjust(hspace=0.275)

ax =fig.add_subplot(211)
ax.plot(dt['Strikes'], dt['Vol_Smile'], label = 'Cubic Spline Smile')
ax.scatter(strike_vec, vol_smile, facecolor = 'orange', linewidth=2.1, label = 'Quoted Smile')
ax.set_title('Vol-Smile for EURNOK')
ax.set_xlabel('Strike')
ax.set_ylabel('Volatility')
ax.legend()
ax.grid(True)


ax = fig.add_subplot(212)
ax.plot(dt['Strikes'], dt['OptionPrice'], label = 'Cubic Spline Prices')
ax.plot(dt['Strikes'], dt['InnerValue'], label = 'Option: Inner Value')
ax.scatter(strike_vec, option_price_vec, color='orange', linewidth=2.1, label='Quoted Prices')
ax.set_title('Option Price for EURNOK Call option along the volatility smile')
ax.set_xlabel('Strike')
ax.set_ylabel('Option Price (NOK/Option)')
ax.legend()
ax.grid(True)

plt.show()

print(dt)




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