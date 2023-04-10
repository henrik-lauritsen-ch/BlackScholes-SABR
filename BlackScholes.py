"""
 FX Vanilla option tools: Garman Kohlhagen + Greeks

"""

import math as m
import Utility as u 
import enum
import unittest


class OptionType(enum.Enum):
    Put = 1
    Call = 2


class GarmanKohlhagen:
    def __init__(self, spot, strike, expiryTerm, depositDomestic, depositForeign, volatility):
        self._spot = spot
        self._strike = strike
        self._expiryTerm = expiryTerm
        self._depositDomestic = depositDomestic
        self._depositForeign = depositForeign
        self._volatility = volatility


class Vanilla(GarmanKohlhagen):  
    def __init__(self, spot, strike, expiryTerm, depositDomestic, depositForeign, volatility):
        super().__init__(spot, strike, expiryTerm, depositDomestic, depositForeign, volatility)

    def Getd1(self, volatility):     
        s0 = self._spot
        K = self._strike   
        r = self._depositDomestic
        d = self._depositForeign
        vol = volatility
        T = self._expiryTerm
        
        if (s0>0 and K>0 and vol>0 and T>0):
            d1 = (m.log(s0/K) + (r - d + 0.5*vol*vol)*T)/(vol*m.sqrt(T))
        else:
            d1 = 0.0       
        return d1
        
    
    def Getd2(self, volatility):
        return self.Getd1(volatility) - volatility*m.sqrt(self._expiryTerm)

    def GetOptionValueSVO(self, strike, volatility, optionType: OptionType):
        
        self._strike = strike
        s0 = self._spot        
        r = self._depositDomestic
        d = self._depositForeign
        T = self._expiryTerm
        d1 = self.Getd1(volatility)
        d2 = self.Getd2(volatility)
        
        sign = 1.0
        if (optionType == OptionType.Put):
            sign = -1.0
        
        return sign * (m.exp(-d * T) * s0 * u.norm().cdf(sign * d1) - m.exp(-r * T) * strike * u.norm().cdf(sign * d2))

    
    def GetBaseOptionValue(self, optionType: OptionType, volatility):        
        return self.GetOptionValueSVO(self._strike, volatility, optionType)


    def GetOptionValue(self, optionType: OptionType):        
        return self.GetBaseOptionValue(optionType, self._volatility)


    def GetOptionValue(self, optionType: OptionType):        
        return self.GetBaseOptionValue(optionType, self._volatility)

        
    def ObjectFuncImpliedVol(self, volatility):                
        return self.GetBaseOptionValue(self._optionType, volatility) - self._targetValue        
        
        
    def GetImpliedVolatility(self, targetValue: float, optionType: OptionType):        
        low = 0.001
        high = 0.90
        self._optionType = optionType
        self._targetValue = targetValue
        
        return u.bisection(self.ObjectFuncImpliedVol, low, high, 0.000000001, 50)


    def GetDomesticSpotDelta(self, optionType: OptionType):        
            
        q = self._depositForeign
        sign = 1.0        
        if (optionType == OptionType.Put):
            sign = -1.0
    
        d1 = self.Getd1(self._volatility)
        signPhi_d1 = u.norm().cdf(sign * d1)

        return sign * m.exp(-q * self._expiryTerm) * signPhi_d1


    def GetGamma(self):
        retval = 0
        d1 = self.Getd1(self._volatility)
        phi_d1 = u.norm().pdf(d1)
        variance = self._volatility * self._volatility * self._expiryTerm
        rootVariance = m.sqrt(variance)

        if ((self._volatility > 0.0) and (self._expiryTerm > 0.0) and (self._spot > 0.0)):
            retval = 1.0 / (rootVariance * self._spot) * m.exp(- self._depositForeign * self._expiryTerm) * phi_d1
        else:
            raise ValueError("Expiry term + volatility + spot needs to be positive - GarmanKohlhagen->Vanilla->GetGamma")

        return retval


    def GetDualDelta(self, optiontype: OptionType):
        
        d2 = self.Getd2(self._volatility)
        
        sign = 1.0
        if (optiontype == OptionType.Put):
            sign = -1.0
        
        return -sign*m.exp(-self._depositDomestic*self._expiryTerm)*u.norm().cdf(sign*d2)


    def GetDualGamma(self):
        d2 = self.Getd2(self._volatility)
        try:
            return m.exp(-self._depositDomestic*self._expiryTerm)*u.norm().pdf(d2)/(self._strike*self._volatility*m.sqrt(self._expiryTerm))
        except:     
            raise ValueError('Divide with zero: GarmanKohlhagen->Vanilla->GetDualGamma')
    
    
    def GetDualVega(self):
        
        # Dual Vega = dVega/dStrike
        d1 = self.Getd1(self._volatility)
        return self._spot/self._strike*m.exp(-self._depositForeign*self._expiryTerm)*d1/self._volatility*u.norm().pdf(d1)
    
    
    def GetVega(self):
        
        d1 = self.Getd1(self._volatility)
        phi_d1 = u.norm().pdf(d1)

        if (self._expiryTerm > 0.0):
            return self._spot * m.exp(- self._depositForeign * self._expiryTerm) * m.sqrt(self._expiryTerm) * phi_d1
        else:
            raise ValueError("Expiry term + volatility needs to be positive - GarmanKohlhagen->Vanilla->GetVega")



    def GetVolga(self):
        retval = 0
        d1 = self.Getd1(self._volatility)
        d2 = self.Getd2(self._volatility)
        phi_d1 = u.norm().pdf(d1)

        if (self._volatility > 0.0 and self._expiryTerm > 0.0):
            retval = self._spot / self._volatility * m.exp(- self._depositForeign * self._expiryTerm) * m.sqrt(self._expiryTerm) * d1 * d2 * phi_d1
        else:
            raise ValueError("Expiry term + volatility needs to be positive - GarmanKohlhagen->Vanilla->GetVolga")

        return retval


    def GetVanna(self):
        retval = 0
        d1 = self.Getd1(self._volatility)
        d2 = self.Getd2(self._volatility)
        phi_d1 = u.norm().pdf(d1)

        if (self._volatility > 0.0 and self._expiryTerm > 0.0):
            retval = -m.exp(-self._depositForeign * self._expiryTerm) * d2 / self._volatility * phi_d1
        else:
            raise ValueError("Expiry term + volatility needs to be positive - GarmanKohlhagen->Vanilla->GetVanna")

        return retval


    def GetTheta(self, optiontype: OptionType):
        
        retval = 0
        if (self._expiryTerm > 0.0):
            
            sign = 1.0            
            if (optiontype == OptionType.Put):
                sign = -1.0

            d1 = self.Getd1(self._volatility)
            d2 = self.Getd2(self._volatility)
            
            Phi_d1 = u.norm().cdf(sign * d1)
            phi_d1 = u.norm().pdf(d1)
            Phi_d2 = u.norm().cdf(sign * d2)

            retval = (-self._spot * m.exp(-self._depositForeign * self._expiryTerm) * phi_d1 * self._volatility / (2 * m.sqrt(self._expiryTerm)) 
            + sign * self._depositForeign * self._spot * Phi_d1 * m.exp(-self._depositForeign * self._expiryTerm) 
            - sign * self._depositDomestic * self._strike * m.exp(-self._depositDomestic * self._expiryTerm) * Phi_d2)

        else:
            raise ValueError("Expiry term + volatility needs to be positive - GarmanKohlhagen->Vanilla->GetTheta")

        return retval



#//     Unit-Test: Garmankohlhagen->Vanilla Class
class TestBSMethods(unittest.TestCase):

    print('Running BlackScholes testing ')
    bs = Vanilla(45.451, 46, 0.876, 0.054, 0.1, 0.18)
    
    def test_GetOptionValue(self):    
        self.assertEqual(round(self.bs.GetOptionValue(OptionType.Put), 6), round(4.12504280568063, 6))

    def test_GetDomesticSpotDelta(self):
        self.assertEqual(round(self.bs.GetDomesticSpotDelta(OptionType.Call), 12), round(0.376084301042287, 12))
        
    def test_GetGamma(self):
        self.assertEqual(round(self.bs.GetGamma(), 8), round(0.0465248841826033, 8))

    def test_GetVega(self):
        self.assertEqual(round(self.bs.GetVega(), 12), round(15.1547507432278, 12))
    
    def test_GetVanna(self):
        self.assertEqual(round(self.bs.GetVanna(), 10), round(0.781155007052361, 10))
        
    def test_GetVolga(self):        
        self.assertEqual(round(self.bs.GetVolga(), 10), round(7.51731525857879, 10))
                
    def test_GetTheta(self):        
        self.assertEqual(round(self.bs.GetTheta(OptionType.Put), 6), round(-2.46333519727929, 6))    
    
    def test_GetImpliedVolatility(self):        
        self.assertEqual(round(self.bs.GetImpliedVolatility(4.206984, OptionType.Put), 4), 0.1854)
        
    def test_GetDualDelta(self):
        self.assertEqual( round(self.bs.GetDualDelta(OptionType.Call), 12), round(-0.330524888341401, 12))
        
    def test_GetDualGamma(self):
        self.assertEqual( round(self.bs.GetDualGamma(), 12), round(0.0454209823850242, 12))

if __name__ == '__main__':
    unittest.main()
    