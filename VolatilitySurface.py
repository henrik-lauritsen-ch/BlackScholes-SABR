"""
 FX Vanilla option tools: Strike From Delta, Forward

"""
import BlackScholes as bs
import math
import Utility as u
import unittest
import numpy as np
import StrikeFromDelta as sfd



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
        sd = sfd.StrikeFromDelta(self._spot, self._domesticDeposit, self._foreignDeposit, self._expiryTerm)
        atm_strike = sd.GetATMStrike(self._ATMVol)
        moneyness_vec = sd.GetLogMoneynessStrikeVector(smile_vec)        
        x = math.log(strike/atm_strike)

        return self._volatilityInterpolation.GetInterpolatedValue(x, moneyness_vec, smile_vec)


class SABRVolSurface(FXVolSurface):
    
    def __init__(self, spot, domesticDeposit, foreignDeposit, expiryTerm, volatilitySmile, beta=0.85):
        super().__init__(spot, domesticDeposit, foreignDeposit, expiryTerm, volatilitySmile)

        self._beta = beta
        sd = sfd.StrikeFromDelta(spot, domesticDeposit, foreignDeposit, expiryTerm)
        self._strikes = sd.GetStrikeVector(volatilitySmile)

        if (expiryTerm <= 1.0):
            self._vovol = 1.0
        elif (expiryTerm <= 3.0):
            self._vovol = 0.5
        else:
            self._vovol = 0.25
        
        self._vovolmin = 0
        self._vovolmax = 0
        
        self._alpha = volatilitySmile[2]        
        self._corr = 0.9999
        self._alphamin = 0
        self._alphamax = 0
        self._sabrrr25 = 0
        
        
    def GetVolatilityFromSmile(self, strike, smile_vec):
                
        # r = self.SabrCalibration(smile_vec)
        # Missing: here we need to map r-vec back to interval
        corr = self._corr #r[0]
        vovol = self._vovol #r[1]
        alpha = self._alpha #r[2]
                        
        return self.SabrImpliedVol(strike, corr, vovol, alpha, self._beta)


    def SabrImpliedVol(self, strike, corr, vovol, alpha, beta):        
        return self.I0_JObloj(strike, corr, vovol, alpha, beta)*(1 + self.I1_Hagan(strike, corr, vovol, alpha, beta)*self._expiryTerm)


    # The I0_JObloj(), the I0() method from "Fine-tune your smile - correction to Hagan et. al." by Jan Oblój
    def I0_JObloj(self, strike, corr, vovol, alpha, beta):               

        forward = sfd.ForwardContinuousDeposit(self._spot, self._domesticDeposit, self._foreignDeposit, self._expiryTerm)
        
        if (forward == strike):
            retval = alpha * math.pow(strike, (beta - 1))
        elif ((vovol == 0.0) and (beta != 1.0)):
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
            
        forward = sfd.ForwardContinuousDeposit(self._spot, self._domesticDeposit, self._foreignDeposit, self._expiryTerm)
        
        return ((1.0 - beta) * (1.0 - beta) / 24.0 * alpha * alpha / math.pow(forward * strike, 1.0 - beta) + 
                1.0 / 4.0 * corr * beta * vovol * alpha / math.pow(forward * strike, (1.0 - beta) / 2.0) + 
                (2.0 - 3.0 * corr * corr) / 24.0 * vovol * vovol)


    # This method established the differencence between the SABR implied vol and the ATM volatility for a given alpha    
    def FirstGuessAlphaMax(self, x):
    
        # Correlationen set to 1.0 in order to solve for the upper bound on alpha (could have have been -1 also)
        # See equation 9 in "Fitting the smile - Smart parameters for SABR and Heston" 
        # by Pierre Gauthier and Pierre-Yves H. Rivaille, 2009 (PP)
    
        # forward = sfd.ForwardContinuousDeposit(self._spot, self._domesticDeposit, self._foreignDeposit, self._expiryTerm)       
        ATMStrike = self._strikes[2]        
        return self.SabrImpliedVol(ATMStrike, 0.9999, self._vovol, x, self._beta) - self._ATMVol
    

    # Input correlation and the difference between the SABR 25 delta risk-reversal and the quoted 25 delta risk
    # reversal is returned. This method is applied for a first guess of the SABR correlation parameter
    def FirstGuessCorrelation(self, x):

        self._sabrrr25 = (self.SabrImpliedVol(self._strikes[3], x, self._vovol, self._alpha, self._beta) 
                    - self.SabrImpliedVol(self._strikes[1], x, self._vovol, self._alpha, self._beta))
        
        return self._sabrrr25 - self._rr25

    def SABRFirstGuess(self) -> None:
        
        # Initial guess and bounds on Max_Alpha given vovol found above and corr = 1.0 we find the upper bound of
        # alpha. We solve equation 9 in PP to establish this bound. Note, that as vovol is a guess then alpha may be
        # larger than the alpha established as the max alpha. For that reason we chooes alphamax = alpha *1.2
        self._alpha = u.bisection(self.FirstGuessAlphaMax, self._ATMVol/2.0, self._ATMVol*3.0, 0.0000001, 25)
        self._alphamin = self._alpha * 0.8       
        self._alphamax = self._alpha * 1.2
               
        self._corr=0
        while (self._corr==0):
        
            try:                
                self._corr = u.bisection(self.FirstGuessCorrelation, -0.9999, 0.99999, 0.00001, 10) 
                print(self._corr)
            except:                
                if (abs(self._sabrrr25) < abs(self._rr25)):
                    self._vovol *= 1.2
                else:
                    self._vovol /= 1.2
                print(self._corr)
        
        self._vovolmin = self._vovol/2.0
        self._vovolmax = self._vovol*2
        
        pass


    def SabrCalibration(self, smile_vec):
        # Missing Build calibration
        return np.array([0, 0.1, 0.1])

vs_vec = np.array([0.14852, 0.10042, 0.0973, 0.11582, 0.23732])
sabr = SABRVolSurface(85, 0.012, 0.00523, 7/365.0, vs_vec)

sabr.SABRFirstGuess()

print(sabr._alpha)
print(sabr._corr)
print(sabr._vovol)


# print(sabr._alphamax)
# print(sabr._alphamin)
# print(sabr._vovolmax)
# print(sabr._vovolmin)
# print(sabr.SabrImpliedVol(85.16357197805945, sabr._corr, sabr._vovol, sabr._alpha, sabr._beta))
# print(sabr.SabrImpliedVol(85.16357197805945, 0.9999, sabr._vovol, sabr._alpha, sabr._beta))
 


#//     Unit-Test: Volatility Surface
class Test_VolSurface(unittest.TestCase):
   
    def test_FXVolSurface_CSI(self):
        vs_vec = np.array([0.117885, 0.1191, 0.1300, 0.1501, 0.174995])        
        csi = u.CubicSplineInterpolation(True)
        vs = FXVolSurface(85.3678, 0.012, 0.0053, 61/365.0, vs_vec, csi)
        self.assertEqual(round(vs.GetVolatility(90.123), 14), round(0.154975045068546, 14))
        self.assertEqual(round(vs.GetVolatility(94.123), 14), round(0.174995, 14))
        
if __name__ == '__main__':
    unittest.main()
