# BlackScholes-SABR
Functionality to value FX options under Garman Kohlhagen and SABR


## BlackScholes.py
**Classes:**
```
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
- InverseCdf(), inverse normal P. J. Acklam.
- Moro(), inverse normal Moro
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
```

```
Methods:

- ForwardContinuousDeposit(spot, domesticDeposit, foreignDeposit, expiryTerm)
```

## SABR.py
**Classes/Methods:**


## VisualizeVolatilitySurfaceFunctionality.py
The purpose of this library is to show application of the different methods implemented for FX Options

 ![Smile](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/bss_fx_smile2_cubic_spline.png)
