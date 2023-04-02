# BlackScholes-SABR
Functionality to value FX options under Garman Kohlhagen and SABR


## BlackScholes.py

This library contains the standard vanilla option pricing formula applied for FX trading (Garman-Kohlhagen) along with standard Greeks.

The Garman-Kohlhagen formula for a vanilla Call:

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/gk_call.png" width=60% height=60%>



Note: For a currency cross CCY1/CCY2 the second currency is considered the domestic currency and the first currency the foreign!

**Classes:**
```
class Option: Option base class

class Vanilla(Option): Garman-Kohlhagen method along with all Greeks

class Exotic(Option): Garman-Kohlhagen version for 1st generation FX exotics

class OptionType(enum.Enum): Put/Call
```

## Utility.py
**Classes/Methods:**

```
class norm: Standard normal density and distribution functions along with inverse normal

class Interpolation: Piecewiese linear and Cubic-Splie 
```

```
Methods:
- RealAxisToIntervalAB()
- IntervalABToRealAxis()
- FindIndex()
- Bisection()
```

## VolatilitySurface.py
**Classes/Methods:**

```
class FXVolSurface: Volatility surface that contains a smile for one Maturity only 
(10dput, 10dput, atm, 10dcall, 10dcall). Constructor takes Interpolation(). Class 
returns volatility given strike.

class SABRVolSurface(FXVolSurface): SABR volatility Surface. Calibrates smile to SABR 
and returns implied SABR-Vol given strike.

```


## StrikeFromDelta.py
**Classes/Methods:**

```
class StrikeFromDelta: Method for calculating Strike-from-delta and ATM-strike
```
```
Methods:
- ForwardContinuousDeposit(spot, domesticDeposit, foreignDeposit, expiryTerm)
```



## The SABR Model
<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_model_v10.png" width=65% height=65%>




## VisualizeVolatilitySurfaceFunctionality.py
The purpose of this library is to show application of the different methods implemented for FX Options

### USD/JPY SABR Calibration:
Implementation: Jan Obloj, Fine-Tune Your Smile Correction to Hagan et al

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/USDJPY_Obloj_1.png" width=60% height=60%>

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/USDJPY_Obloj_2.png" width=60% height=60%>



### ZAR/JPY SABR Calibration:
Implementation: Jan Obloj, Fine-Tune Your Smile Correction to Hagan et al

![Smile](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_ZARJPY.png)

### ZAR/JPY SABR Calibration:
Implementation: Jan Obloj + Benaim/Dodgson/Kainth Smile Extrapolation

![Smile](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_ZARJPY_SE.png)



### CHF/DKK SABR Calibration:
Implementation: Jan Obloj + Benaim/Dodgson/Kainth Smile Extrapolation

![Smile](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_CHFDKK.png)

