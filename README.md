# BlackScholes-SABR
Functionality to value FX options under Garman Kohlhagen and SABR


## BlackScholes.py

This library contains the standard vanilla option pricing formula applied for FX trading, Garman-Kohlhagen, along with standard Greeks.

The formula Garman-Kohlhagen formula for a vanilla Call option:

![Garman-Kohlhagen](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/garman-kohlhagen.png)


**Classes:**
```

.... Garman Kohlhagen ..... EURUSD 1. + 2. ccy was heisst long/short?

class Vanilla:

- GetOptionValue(optionType)
- GetBaseOptionValue(optionType: OptionType, volatility)
- Getd1(volatility)
- Getd2(volatility)
- GetImpliedVolatility(targetValue, optionType)
- GetDomesticSpotDelta(optionType)
- GetGamma()
- GetVega()
- GetVolga()
- GetVanna()
- GetTheta()
```

```
class OptionType(enum.Enum):

- Put
- Call
```
## Utility.py
**Classes/Methods:**


```
class norm:

- pdf(), standard normal probability density function
- cdf(), standard normal cumulative distribution function
- cdfM(), standard normal cumulative distribution function (method M)
- cdfI(), standard normal cumulative distribution function (method I)
- InverseCdf(), inverse normal, P. J. Acklam.
- Moro(), inverse normal, Moro
```

```
class Interpolation:

- LinearInterpolation()
- PiecewiseLinearInterpolation()
- CubicSplineInterpolation()
```

```
Methods:

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

