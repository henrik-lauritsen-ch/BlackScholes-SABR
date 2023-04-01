# BlackScholes-SABR
Functionality to value FX options under Garman Kohlhagen and SABR


## BlackScholes.py

This library contains the standard vanilla option pricing formula applied for FX trading, Garman-Kohlhagen along with standard Greeks.

The formula Garman-Kohlhagen formula for a vanilla Call option:

![Garman-Kohlhagen](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/gk_call.png)

with

![Garman-Kohlhagen](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/gk_d1d2.png)

and where 

![Garman-Kohlhagen](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/gk_discountfactor.png)

![Garman-Kohlhagen](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/gk_forward.png)

S: FX Spot

K: Strike

t: today (can be set to zero)

T: Maturity term (in years)

r: Domestic deposit (continous, second currency)

q: Foreign deposit (continous, first currency)

N(): Standard cumulative normal distribution function

sigma: traded implied volatility

**Classes:**
```
class Vanilla: Garman-Kohlhagen method along with all Greeks
```

```
class OptionType(enum.Enum): Put/Call
```

## Utility.py
**Classes/Methods:**


```
class norm: Standard normal dentisity and distribution functions along with inverse normal
```

```
class Interpolation: Piecewiese linear and Cubic-Splie 
```

```
Methods:
- RealAxisToIntervalAB()
- IntervalABToRealAxis
- FindIndex()
- Bisection()
```

## VolatilitySurface.py
**Classes/Methods:**

```
class StrikeFromDelta:

- GetATMStrike(self, volatility)
- GetStrikeFromDomesticDelta(self, delta, optiontype: bs.OptionType, volatility)
- GetStrikeVector(self, volSmile)
- GetLogMoneynessStrikeVector(self, volSmile)
```

```
Methods:

- ForwardContinuousDeposit(spot, domesticDeposit, foreignDeposit, expiryTerm)
```

## SABR.py
**Implementation and Calibration of the SABR model:**

 ![Smile](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr2_equations.png)

**Classes/Methods:**

## VisualizeVolatilitySurfaceFunctionality.py
The purpose of this library is to show application of the different methods implemented for FX Options

### USD/JPY SABR Calibration:
![Smile](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_USDJPY.png)

### ZAR/JPY SABR Calibration:
![Smile](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_ZARJPY.png)

### CHF/DKK SABR Calibration:
![Smile](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_CHFDKK.png)

