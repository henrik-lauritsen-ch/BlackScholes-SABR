"""
 FX Vanilla option tools: Strike From Delta, Forward

"""
import math as m
import Utility as u
import unittest
import numpy as np
import StrikeFromDelta as sfd
import scipy.optimize as so


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
        self._sd = sfd.StrikeFromDelta(spot, domesticDeposit, foreignDeposit, expiryTerm)
        self._strikes = np.nan
    
    
    def SetVolatilitySmile(self, volatilitysmile) -> None:    
        self.UpdateVolatilities(volatilitysmile)
        self.CalcStrikeVector()
        pass
    
    
    def UpdateVolatilities(self, volatilitysmile) -> None:
        self._volatilitySmile = volatilitysmile
        self._ATMVol = volatilitysmile[2]
        self._rr25 = volatilitysmile[3] - volatilitysmile[1]
        pass
    
    
    def GetVolatility(self, strike):        
        return self.GetVolatilityFromSmile(strike, self._volatilitySmile)


    def GetVolatilityFromSmile(self, strike, smile_vec):
        
        self.SetVolatilitySmile(smile_vec)
        atm_strike = self._strikes[2]
        moneyness_vec = self._sd.GetLogMoneynessStrVec(self._strikes)
        x = m.log(strike/atm_strike)

        return self._volatilityInterpolation.GetInterpolatedValue(x, moneyness_vec, smile_vec)


    def CalcStrikeVector(self) -> None:        
        self._strikes = self._sd.GetStrikeVector(self._volatilitySmile)        
        pass


class SABRVolSurface(FXVolSurface):
    
    def __init__(self, spot, domesticDeposit, foreignDeposit, expiryTerm, volatilitySmile, beta=0.85):
        super().__init__(spot, domesticDeposit, foreignDeposit, expiryTerm, volatilitySmile)

        self._beta = beta              
        self._alpha0 = np.nan
        self._alpha = np.nan
        self._alphamin = np.nan
        self._alphamax = np.nan
        self._corr0 = np.nan
        self._corr = np.nan
        self._vovol0 = np.nan
        self._vovol = np.nan        
        self._vovolmin = np.nan
        self._vovolmax = np.nan
        
        self._sabrrr25 = np.nan        
        
        self.CalcStrikeVector()
        self._calibrationWeights = np.array([1.0, 1.0, 2.0, 1.0, 1.0])              
        self.SabrCalibration()
        
        
    def GetVolatilityFromSmile(self, strike, smile_vec):
        
        if (smile_vec != self._volatilitySmile).any():        
            self.SetVolatilitySmile(smile_vec)        
            self.SabrCalibration()
            
        return self.SabrImpliedVol(strike, self._alpha, self._corr, self._vovol, self._beta)
  

    def SabrImplVolFwd(self, strike, forward, alpha, corr, vovol, beta):
        return self.I0_JObloj(strike, forward, alpha, corr, vovol, beta)*(1.0 + self.I1_Hagan(strike, forward, alpha, corr, vovol, beta)*self._expiryTerm)


    def SabrImpliedVol(self, strike, alpha, corr, vovol, beta):  
        forward = sfd.ForwardContinuousDeposit(self._spot, self._domesticDeposit, self._foreignDeposit, self._expiryTerm)      
        return self.SabrImplVolFwd(strike, forward, alpha, corr, vovol, beta)


    # The I0_JObloj(), the I0() method from "Fine-tune your smile - correction to Hagan et. al." by Jan Oblój
    def I0_JObloj(self, strike, forward, alpha, corr, vovol, beta):
        
        x = self.GetI0x(strike, forward)        
        
        if (strike == forward):
            retval = alpha * pow(strike, (beta - 1.0))
        elif ((vovol == 0.0) and (beta != 1.0)):
            retval = x * alpha * (1.0 - beta) / (pow(forward, (1.0 - beta)) - pow(strike, (1.0 - beta)))
        else:            
            z = self.GetI0z(strike, forward, alpha, vovol, beta)            
            denominator = m.log((m.sqrt(1.0 - 2.0 * corr * z + z * z) + z - corr) / (1.0 - corr))
            retval = vovol * x / denominator

        return retval


    def GetI0z(self, strike, forward, alpha, vovol, beta):
        
        x = self.GetI0x(strike, forward)
        
        if (beta == 1.0):
            z = vovol * x / alpha
        else:
            z = (vovol / alpha) * (pow(forward, (1.0 - beta)) - pow(strike, (1.0 - beta)))/(1.0 - beta)

        return z
    
    
    def GetI0x(self, strike, forward):
        
        try:
            return m.log(forward/strike)
        except:        
            raise ValueError('Divide with zero or ln(0): SABRVolSur->GEtI0x')
        
        
    # The I1_Hagan(), second term of the implied volatility formula (eq. 2.17a) in "Managing Smile Risk" by Hagan et. al
    def I1_Hagan(self, strike, forward, alpha, corr, vovol, beta):
        
        return (pow((beta - 1.0), 2.0)/24.0 * pow(alpha, 2.0) / pow(forward * strike, 1.0 - beta) + 
                1.0/4.0 * corr * vovol * alpha * beta / pow(forward * strike, (1.0 - beta)/2.0) + 
                (2.0 - 3.0 * corr * corr) / 24.0 * vovol * vovol)


    # This method established the differencence between the SABR implied vol and the ATM volatility for a given alpha    
    def FirstGuessAlphaMax(self, x):
    
        # Correlationen set to 1.0 in order to solve for the upper bound on alpha (could have have been -1 also)
        # See equation 9 in "Fitting the smile - Smart parameters for SABR and Heston" 
        # by Pierre Gauthier and Pierre-Yves H. Rivaille, 2009 (PP)        
        ATMStrike = self._strikes[2]        
        return self.SabrImpliedVol(ATMStrike, x, 0.9999, self._vovol0, self._beta) - self._ATMVol
    

    # Input correlation and the difference between the SABR 25 delta risk-reversal and the quoted 25 delta risk
    # reversal is returned. This method is applied for a first guess of the SABR correlation parameter
    def FirstGuessCorrelation(self, x):

        self._sabrrr25 = (self.SabrImpliedVol(self._strikes[3], self._alpha0, x, self._vovol0, self._beta) 
                         -self.SabrImpliedVol(self._strikes[1], self._alpha0, x, self._vovol0, self._beta))
        
        return self._sabrrr25 - self._rr25


    def SABRFirstGuess(self) -> None:
                
        if (self._expiryTerm <= 14.0/365.0):
            self._vovol0 = 5.0
        elif (self._expiryTerm <= 31.0/365.0):
            self._vovol0 = 2.0
        elif (self._expiryTerm <= 1.0):
            self._vovol0 = 1.0
        elif (self._expiryTerm <= 3.0):
            self._vovol0 = 0.5
        else:
            self._vovol0 = 0.25
                
                
        # Initial guess and bounds on Max_Alpha given vovol found above and corr = 1.0 we find the upper bound of
        # alpha. We solve equation 9 in PP to establish this bound. Note, that as vovol is a guess then alpha may be
        # larger than the alpha established as the max alpha. For that reason we chooes alphamax = alpha *1.2
        self._alpha0 = u.bisection(self.FirstGuessAlphaMax, self._ATMVol/2.0, self._ATMVol*3.0, 0.00001, 10)
        self._alphamin = self._alpha0 * 0.7       
        self._alphamax = self._alpha0 * 1.2
               
        self._corr0=0
        while (self._corr0==0.0):
        
            try:                
                self._corr0 = u.bisection(self.FirstGuessCorrelation, -0.9999, 0.99999, 0.00001, 10)                 
            except:                
                if (abs(self._sabrrr25) < abs(self._rr25)):
                    self._vovol0 *= 1.2
                else:
                    self._vovol0 /= 1.2
        
        self._vovolmin = self._vovol0/2.0
        self._vovolmax = self._vovol0*3.0
        
        pass


    def SABRCalibObject(self, v):
        
        s = 0.0
        for i in range(5):
            
            sabrvol = self.SabrImpliedVol(self._strikes[i], u.RealAxisToIntervalAB(v[0], self._alphamin, self._alphamax),
                                                            u.RealAxisToIntervalAB(v[1], -0.9999, 0.9999),
                                                            u.RealAxisToIntervalAB(v[2], self._vovolmin, self._vovolmax),
                                                            self._beta)
            s += pow(sabrvol - self._volatilitySmile[i], 2) * self._calibrationWeights[i]
            
        return s
        

    def SabrCalibration(self) -> None:
        
        # Pre-calibation - First guess on a solution
        self.SABRFirstGuess()
        
        # Calibrate, parameters mapped to R3
        x0 = np.array([u.IntervalABToRealAxis(self._alpha0, self._alphamin, self._alphamax),
                       u.IntervalABToRealAxis(self._corr0, -0.9999, 0.9999), 
                       u.IntervalABToRealAxis(self._vovol0, self._vovolmin, self._vovolmax)])
        x_min = so.minimize(self.SABRCalibObject, x0, method='Powell', tol=0.0000001)
        # print(x_min)
        
        # Map result back tp min/max intervals
        self._alpha = u.RealAxisToIntervalAB(x_min.x[0], self._alphamin, self._alphamax)
        self._corr = u.RealAxisToIntervalAB(x_min.x[1], -0.9999, 0.9999)
        self._vovol= u.RealAxisToIntervalAB(x_min.x[2], self._vovolmin, self._vovolmax)
        
        pass


    
       
        
#//     Unit-Test: Volatility Surface
class Test_VolSurface(unittest.TestCase):
       
    def test_FXVolSurface_CSI(self):
        vs_vec = np.array([0.117885, 0.1191, 0.1300, 0.1501, 0.174995])        
        csi = u.CubicSplineInterpolation(True)
        vs = FXVolSurface(85.3678, 0.012, 0.0053, 61/365.0, vs_vec, csi)
        self.assertEqual(round(vs.GetVolatility(90.123), 14), round(0.154975045068546, 14))
        self.assertEqual(round(vs.GetVolatility(94.123), 14), round(0.174995, 14))
        
    def test_SABR(self):
        # Todo Freitag!!
        pass 
    
if __name__ == '__main__':
    unittest.main()
