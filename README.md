# BlackScholes-SABR
Functionality to value FX options under Garman Kohlhagen and SABR


## BlackScholes.py
**Classes:**
```
class Vanilla:
```
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
class OptionType(enum.Enum):
```
- Put
- Call

## Utility.py
**Classes/Methods:**
   

## VolatilitySurface.py
**Classes/Methods:**


## SABR.py
**Classes/Methods:**


## VisualizeVolatilitySurfaceFunctionality.py
The purpose of this library is to show application of the different methods implemented for FX Options

 ![Smile](https://github.com/henrik-lauritsen-ch/Pictures/blob/main/bss_fx_smile_cubic_spline.png)
