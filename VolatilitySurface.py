"""
 FX Vanilla option tools: Strike From Delta, Forward

"""
import BlackScholes as bs
import math
import Utility as u
import unittest
import numpy as np


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
        
        
if __name__ == '__main__':
    unittest.main()
