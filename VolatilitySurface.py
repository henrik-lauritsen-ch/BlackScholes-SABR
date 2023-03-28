"""
 FX Vanilla option tools: Strike From Delta, Forward

"""
import BlackScholes as bs
import math
import Utility as u
import unittest
import numpy as np


class FXVolSurface:
    
    def __init__(self, spot, domesticDeposit, foreignDeposit, expiryTerm, volatilitySmile, volatilityInterpolation: u.Interpolation = u.CubicSplineInterpolation(False)):

        self._spot = spot
        self._expiryTerm = expiryTerm
        self._volatilitySmile = volatilitySmile
        self._domesticDeposit = domesticDeposit
        self._foreignDeposit = foreignDeposit
        self._volatilityInterpolation = volatilityInterpolation        
        self._ATMVol = volatilitySmile[2]
        self._rr25 = volatilitySmile[3] - volatilitySmile[1]

    def GetVolatility(self, strike):        
        return self.GetVolatilityFromSmile(strike, self._volatilitySmile)


    def GetVolatilityFromSmile(self, strike, smile_vec):
        sfd = StrikeFromDelta(self._spot, self._domesticDeposit, self._foreignDeposit, self._expiryTerm)
        atm_strike = sfd.GetATMStrike(self._ATMVol)
        moneyness_vec = sfd.GetLogMoneynessStrikeVector(smile_vec)        
        x = math.log(strike/atm_strike)

        return self._volatilityInterpolation.GetInterpolatedValue(x, moneyness_vec, smile_vec)


class SABRVolSurface(FXVolSurface):
    
    def __init__(self, spot, domesticDeposit, foreignDeposit, expiryTerm, volatilitySmile, beta=0.85):
        super().__init__(spot, domesticDeposit, foreignDeposit, expiryTerm, volatilitySmile)

        self._beta = beta
        sfd = StrikeFromDelta(spot, domesticDeposit, foreignDeposit, expiryTerm)
        self._strikes = sfd.GetStrikeVector(volatilitySmile)

        if (expiryTerm <= 1.0):
            self._vovol = 1.0
        elif (expiryTerm <= 3.0):
            self._vovol = 0.5
        else:
            self._vovol = 0.25


    def GetVolatilityFromSmile(self, strike, expiryterm, smile_vec):
                
        r = self.SabrCalibration(smile_vec)
        # Missing: here we need to map r-vec back to interval
        corr = r[0]
        vovol = r[1]
        alpha = r[2]
                        
        return self.SabrImpliedVol(strike, corr, vovol, alpha, self._beta)


    def SabrImpliedVol(self, strike, corr, vovol, alpha, beta):        
        return self.I0_JObloj(strike, corr, vovol, alpha, beta)*(1 + self.I1_Hagan(strike, corr, vovol, alpha, beta)*self._expiryTerm)


    # The I0_JObloj(), the I0() method from "Fine-tune your smile - correction to Hagan et. al." by Jan Oblój
    def I0_JObloj(self, strike, corr, vovol, alpha, beta):               

        forward = ForwardContinuousDeposit(self._spot, self._domesticDeposit, self._foreignDeposit, self._expiryTerm)
        
        if (forward == strike):
            retval = alpha * math.Pow(strike, (beta - 1))
        elif ((vovol == 0.0) and (beta!=1.0)):
            retval = math.log(forward / strike) * alpha * (1 - beta) / (math.pow(forward, (1 - beta)) - math.pow(strike, (1 - beta)))
        else:
                  
            if (beta == 1.0):
                z = vovol * math.log(forward / strike) / alpha
            else:
                z = vovol / alpha * 1 / (1 - beta) * (math.pow(forward, (1.0 - beta)) - math.pow(strike, (1.0 - beta)))

            denominator = math.log((z - corr + math.sqrt(1 - 2 * corr * z + z * z)) / (1 - corr))
            retval = vovol * math.log(forward / strike) / denominator

        return retval


    # The I1_Hagan(), second term of the implied volatility formula (eq. 2.17a) in "Managing Smile Risk" by Hagan et. al
    def I1_Hagan(self, strike, corr, vovol, alpha, beta):
            
        forward = ForwardContinuousDeposit(self._spot, self._domesticDeposit, self._foreignDeposit, self._expiryTerm)
        
        return ((1.0 - beta) * (1.0 - beta) / 24.0 * alpha * alpha / math.pow(forward * strike, 1.0 - beta) + 
                1.0 / 4.0 * corr * beta * vovol * alpha / math.pow(forward * strike, (1.0 - beta) / 2.0) + 
                (2.0 - 3.0 * corr * corr) / 24.0 * vovol * vovol)


    # This method established the differencence between the SABR implied vol and the ATM volatility for a given alpha    
    def FirstGuessAlphaMax(self, x):
    
        # Correlationen set to 1.0 in order to solve for the upper bound on alpha (could have have been -1 also)
        # See equation 9 in "Fitting the smile - Smart parameters for SABR and Heston" 
        # by Pierre Gauthier and Pierre-Yves H. Rivaille, 2009 (PP)
    
        forward = ForwardContinuousDeposit(self._spot, self._domesticDeposit, self._foreignDeposit, self._expiryTerm)
        return self.SabrImpliedVol(forward, self._expiryTerm, 1.0, self._vovol, x, self._beta) - self._ATMVol
    

    # Input correlation and the difference between the SABR 25 delta risk-reversal and the quoted 25 delta risk
    # reversal is returned. This method is applied for a first guess of the SABR correlation parameter
    def FirstGuessCorrelation(self, x):

        _sabrrr25 = (self.SabrImpliedVol(self._strikes[3], self._expiryTerm, x, self._vovol, _alpha, self._beta) 
                    - self.SabrImpliedVol(self._strikes[1], self._expiryTerm, x, self._vovol, _alpha, self._beta))
        return _sabrrr25 - self._rr25



    def SabrCalibration(self, smile_vec):
        # Missing Build calibration
        return np.array([0, 0.1, 0.1])





def ForwardContinuousDeposit(spot, domesticDeposit, foreignDeposit, expiryTerm):
    return spot * math.exp((domesticDeposit - foreignDeposit) * expiryTerm)


class StrikeFromDelta:
    
    def __init__(self, spot, domesticDeposit, foreignDeposit, expiryTerm):
        self._spot = spot
        self._domesticDeposit = domesticDeposit
        self._foreignDeposit = foreignDeposit
        self._expiryTerm = expiryTerm

    
    def GetATMStrike(self, volatility):
        fwd = ForwardContinuousDeposit(self._spot, self._domesticDeposit, self._foreignDeposit, self._expiryTerm)                
        
        # Below in case Premium-Type should be implemented
        # if (_premiumtype == PremiumType.DOMESTIC):
        #     retval = fwd * math.exp(0.5 * volatility * volatility * self._expiryTerm)
        # elif (_premiumtype == PremiumType.FOREIGN):
        #     retval = fwd * math.exp(-0.5 * volatility * volatility * self._expiryTerm)
        # else:
        #     raise ValueError('GetATMStrike Exception')        
        return fwd * math.exp(0.5 * volatility * volatility * self._expiryTerm)

    
    def GetStrikeFromDomesticDelta(self, delta, optiontype: bs.OptionType, volatility):
        
        sign = -1.0
        if (optiontype == bs.OptionType.Call):        
            sign = 1.0

        delta = abs(delta)

        if (delta >= 1.0):
            delta = 0.999
        
        # Below in case Delta-Type should be implemented
        #if (deltaType = DeltaType:SPOT):
        z = math.exp(self._foreignDeposit*self._expiryTerm)*delta    
        # else:
        #     z = delta
                        
        if (z>=1.0):
            raise ValueError('No solution for this delta and/or these parameters')
        
        norm_inverse = u.norm().InverseCdf(z)
        forward = ForwardContinuousDeposit(self._spot, self._domesticDeposit, self._foreignDeposit, self._expiryTerm)
        
        if (volatility>2):
            raise ValueError('The volatility should be below 200%')
        
        return forward*math.exp(-sign*norm_inverse*volatility*math.sqrt(self._expiryTerm) + 0.5*volatility*volatility*self._expiryTerm)


    def GetStrikeVector(self, volSmile):

        retval = np.zeros(5)       
        retval[0] = self.GetStrikeFromDomesticDelta(0.1, bs.OptionType.Put, volSmile[0])
        retval[1] = self.GetStrikeFromDomesticDelta(0.25, bs.OptionType.Put, volSmile[1])
        retval[2] = self.GetATMStrike(volSmile[2])
        retval[3] = self.GetStrikeFromDomesticDelta(0.25, bs.OptionType.Call, volSmile[3])
        retval[4] = self.GetStrikeFromDomesticDelta(0.1, bs.OptionType.Call, volSmile[4])

        return retval
 

    def GetLogMoneynessStrikeVector(self, volSmile):
        
        strike_vec = self.GetStrikeVector(volSmile)
        atm_strike = self.GetATMStrike(volSmile[2])
        moneyness_vec = np.zeros(5)
        
        for i in range(5):
            moneyness_vec[i] = math.log(strike_vec[i]/atm_strike)
        
        return moneyness_vec



#//     Unit-Test: Volatility Surface
class Test_VolSurface(unittest.TestCase):
    
    sfd = StrikeFromDelta(100, 0.01, 0.02, 1.0)
       
    def test_GetStrikeFromDomDelta(self):
        self.assertEqual(round(self.sfd.GetStrikeFromDomesticDelta(0.23, bs.OptionType.Call, 0.12), 12), round(108.766767822761, 12))

    def test_GetAtmVolatility(self):
        self.assertEqual(round(self.sfd.GetATMStrike(0.11), 10), round(99.6057790988489, 10))

    def test_GetForward(self):
        self.assertEqual(round(ForwardContinuousDeposit(100, 0.01, 0.05, 0.7), 12), round(97.2388366801247, 12))

    def test_GetStrikeVec(self):
        strikes = np.array([9.796265875871027, 10.067098505250692, 10.356101824110898, 10.687697656702378, 11.069582777590423])
        vols = np.array([0.09852, 0.09542, 0.0973, 0.10582, 0.11732])
        sfd2 = StrikeFromDelta(10.3719, 0.00565, 0.01822, 0.194520547945205)
        vcs_strikes = sfd2.GetStrikeVector(vols)
        self.assertEqual(round((strikes - vcs_strikes).sum(), 8), 0.0)

    def test_FXVolSurface_CSI(self):
        vs_vec = np.array([0.117885, 0.1191, 0.1300, 0.1501, 0.174995])        
        csi = u.CubicSplineInterpolation(True)
        vs = FXVolSurface(85.3678, 0.012, 0.0053, 61/365.0, vs_vec, csi)
        self.assertEqual(round(vs.GetVolatility(90.123), 14), round(0.154975045068546, 14))
        self.assertEqual(round(vs.GetVolatility(94.123), 14), round(0.174995, 14))
        
if __name__ == '__main__':
    unittest.main()
