# BlackScholes-SABR
Functionality to value FX options under Garman Kohlhagen and SABR


## BlackScholes.py

This library contains the standard vanilla option pricing formula applied for FX trading (Garman-Kohlhagen) along with standard Greeks.

The Garman-Kohlhagen formula for a vanilla Call:

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/gk_call.png" width=60% height=60%>



Note: For a currency cross CCY1/CCY2 the second currency is considered the domestic currency and the first currency the foreign!

**Classes:**
```
class GarmanKohlhagen: Option base class

class Vanilla(GarmanKohlhagen): Garman-Kohlhagen method along with all Greeks

class OptionType(enum.Enum): Put/Call
```

## ExoticFX.py

**Classes:**
```
class Exotic(GarmanKohlhagen): Garman-Kohlhagen version of 1st generation FX exotics

class ExoticType(enum.Enum): ...
```



## Utility.py
**Classes/Methods:**

```
class norm: Standard normal density and distribution functions along with inverse normal

class Interpolation: Piecewiese linear and Cubic-Spline 
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

class SABRWingVolSurface(SABRVolSurface): SABR volatility Surface. Calibrates smile to SABR.
Extrapolate below 25dPut and above 25dCall using polynomial in prices. The GetVolatility
methods returns implied SABR-Vol given strike.

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

**SABR Dynamics:**

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_model_v10.png" width=65% height=65%>


**SABR Implementation:**

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_implementation_v2.png" width=65% height=65%>

## The SABRWing.py
<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_put_wing1.png" width=65% height=65%>

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_put_wing23.png" width=65% height=65%>

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_put_wing3.png" width=65% height=65%>

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_put_wing4.png" width=65% height=65%>

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_call_wing_extrapolation.png" width=65% height=65%>


Benaim/Dodgson/Kainth suggest ways of dealing with possible arbitrage on the smile and even an alternative solution. These topics have not been taken into account in this Python implementation.



## VisualizeVolatilitySurfaceFunctionality.py
The purpose of this library is to show application of the different methods implemented for FX Options

### USD/JPY SABR Calibration:
Implementation: Jan Obloj, Fine-Tune Your Smile Correction to Hagan et al

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/USDJPY_Obloj_1.png" width=100% height=100%>

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/USDJPY_Obloj_2.png" width=100% height=100%>



### ZAR/JPY SABR Calibration:
Implementation: Jan Obloj, Fine-Tune Your Smile Correction to Hagan et al

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/ZARJPY_Obloj_1.png" width=100% height=100%>

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/ZARJPY_Obloj_2.png" width=100% height=100%>


### ZAR/JPY SABRWing Calibration:
Below we compare the Jan Obloj to the Benaim/Dodgson/Kainth (Wing extrapolation). The implied reisk neutral distribution have caclulated for all four interpolation methods. Obvious that piecewise linear interpolation leads to an unacceptable distribution function.

Clearly a smile like the 2008 ZAR/JPY show the need for more complex volatility modelling.

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_ZARJPY_Wing.png" width=70% height=70%>

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_ZARJPY_PDF.png" width=70% height=70%>





### CHF/DKK SABRWing Calibration:
Implementation: Jan Obloj + Benaim/Dodgson/Kainth (Wing extrapolation)

The need for a more complex volatility model is less imminent for a smile like CHF/DKK as the graphs below show.


<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_CHFDKK_Wing.png" width=100% height=100%>

<img src="https://github.com/henrik-lauritsen-ch/Pictures/blob/main/sabr_CHFDKK_PDF.png" width=100% height=100%>

 
