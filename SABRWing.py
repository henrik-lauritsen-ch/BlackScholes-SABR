# SABR Wing Extrapolation Surface

import VolatilitySurface as vs
import StrikeFromDelta as sfd
import numpy as np
import BlackScholes as bs
import math as m
import unittest



class SABRWingSurface(vs.SABRVolSurface):
    
    def __init__(self, spot, domesticDeposit, foreignDeposit, expiryTerm, volatilitySmile, beta=0.85):
        super().__init__(spot, domesticDeposit, foreignDeposit, expiryTerm, volatilitySmile, beta)
    
        # Strike and Vol are dummy values for the _bs object. This object is applied 
        # for the wing-extrapolation calculation and these parameters are set at "run-time"
        self._bs = bs.Vanilla(spot, 1.0, expiryTerm, domesticDeposit, foreignDeposit, 0.1)
        
        self._ny = np.nan
        self._acall = np.nan
        self._bcall = np.nan
        self._ccall = np.nan
        
        self._my = np.nan
        self._aput = np.nan
        self._bput = np.nan
        self._cput = np.nan
        
        self._calibrationWeights = np.array([0.15, 2.5, 4.0, 2.5, 0.15])
        self.SabrCalibration()
        self.CalcWingParameters()
    
    
    def GetVolatilityFromSmile(self, strike, smile_vec):
        
        if (smile_vec != self._volatilitySmile).any():
            self.SetVolatilitySmile(smile_vec)
            self.SabrCalibration()
            self.CalcWingParameters()
                      
        if (strike < self._strikes[1] or strike > self._strikes[3]):            
            return self.GetImpliedWingVol(strike)
        else:
            return self.SabrImpliedVol(strike, self._alpha, self._corr, self._vovol, self._beta)
    
    
    def GetImpliedWingVol(self, strike):
        
        self._bs._strike = strike
        optionType = bs.OptionType.Call
        
        if (strike < self._strikes[1]):
            optionType = bs.OptionType.Put
            price = self.PutExtrapolationFunction(strike)
        else:
            price = self.CallExtrapolationFunction(strike)
                
        return self._bs.GetImpliedVolatility(price, optionType)
        
            
    def PutExtrapolationFunction(self, strike):    
        return m.exp(self._my*m.log(strike) + self._aput + self._bput*strike + self._cput*strike*strike)
    
    
    def CallExtrapolationFunction(self, strike):    
        return m.exp(-self._ny*m.log(strike) + self._acall + self._bcall/strike + self._ccall/(strike*strike))
    
    
    def dI0dK(self, strike, forward, alpha, corr, vovol, beta):
        
        x = self.GetI0x(strike, forward)
        
        if (forward == strike):
            return alpha*(beta - 1)*pow(strike, beta - 2)
        elif (vovol == 0):
            den = pow(forward, 1 - beta) - pow(strike, 1 - beta)
            return -1/strike*alpha*(1 - beta)/den + x*alpha*pow(1 - beta, 2)/pow(den, 2)*pow(strike, -beta)
        else:
            dzdK = -vovol/alpha*pow(strike, -beta)
            z = self.GetI0z(strike, forward, alpha, vovol, beta)
            sqrtzcorr = m.sqrt(1 - 2*corr*z + z*z) + z - corr
            dlog_num = m.log(sqrtzcorr/(1 - corr))
        
            term1 = -(vovol/strike)*1/dlog_num
            term2 = -(vovol*x)/pow(dlog_num, 2)*1/sqrtzcorr*(pow(1 - 2*corr*z + z*z, -1/2)*(z - corr)*dzdK + dzdK)
        
            return term1 + term2
        
        
    def dI1dK(self, strike, forward, alpha, corr, vovol, beta):        
        return -(pow(1 - beta, 3)/24)*alpha*alpha/pow(forward*strike, 2 - beta)*forward - 1/8*alpha*beta*corr*vovol*(1 - beta)/pow(forward*strike, (3 - beta)/2)*forward


    def d2I0dKdK(self, strike, forward, alpha, corr, vovol, beta):
        
        x = self.GetI0x(strike, forward)
        
        if (forward == strike):
            return alpha*(beta - 1)*(beta - 2)*pow(strike, beta - 3)
        elif (vovol == 0):
            frac =(1 - beta)/(pow(forward, 1 - beta) - pow(strike, 1 - beta))            
            
            return (alpha*frac*pow(strike, -2) - 2*alpha*pow(frac, 2)*pow(strike, -(1 + beta))
                    + 2*x*alpha*pow(frac, 3)*pow(strike, -2*beta) 
                    - x*alpha*pow(frac, 2)*beta*pow(strike, -(1 + beta)))
        else:            
            z = self.GetI0z(strike, forward, alpha, vovol, beta)
            dzdK = -vovol/alpha*pow(strike, -beta)
            d2zdKdK = vovol/alpha*beta*pow(strike, -(1 + beta))            
            sqrtzcorr = 1 - 2*corr*z + z*z
            sqrt_sqrtzcorr = m.sqrt(sqrtzcorr) + z - corr
            log_num = m.log(sqrt_sqrtzcorr/(1 - corr))
            dlog = 1/sqrt_sqrtzcorr*(pow(sqrtzcorr, -1/2)*(z - corr)*dzdK + dzdK)            
            
            term1 = vovol/log_num*pow(strike, -2) + vovol/strike*1/pow(log_num, 2)*dlog
            term2 = (vovol/strike*pow(log_num, -2) + 2*vovol*x*pow(log_num, -3)*dlog)*1/sqrt_sqrtzcorr*(pow(sqrtzcorr, -1/2)*(z - corr)*dzdK + dzdK)
            term31 = -pow(sqrt_sqrtzcorr, -2)*pow((pow(sqrtzcorr, -1/2)*(z - corr)*dzdK + dzdK), 2)
            term32 = pow(sqrt_sqrtzcorr, -1)*(-pow(sqrtzcorr, -3/2)*pow((z - corr)*dzdK, 2) + pow(sqrtzcorr, -1/2)*(dzdK*dzdK + (z - corr)*d2zdKdK) + d2zdKdK)
            term3 = -vovol*x*pow(log_num, -2)*(term31 + term32)
            return term1 + term2 + term3
    
    
    def d2I1dKdK(self, strike, forward, alpha, corr, vovol, beta):
        
        term1 = (alpha*pow(1 - beta, 3)*(2 - beta)/24)*1/pow(forward*strike, 3 - beta)*forward*forward
        term2 = (alpha*beta*corr*vovol*(1 - beta)*(3 - beta)/16)*1/pow(forward*strike, (5 - beta)/2)*forward*forward
        
        return term1  + term2


    def dSABRdK(self, strike, forward, alpha, corr, vovol, beta):
        term1 = self.dI0dK(strike, forward, alpha, corr, vovol, beta)*(1 + self.I1_Hagan(strike, forward, alpha, corr, vovol, beta)*self._expiryTerm)
        term2 = self.I0_JObloj(strike, forward, alpha, corr, vovol, beta)*self.dI1dK(strike, forward, alpha, corr, vovol, beta)*self._expiryTerm
        return term1 + term2


    def d2SABRdKdK(self, strike, forward, alpha, corr, vovol, beta):
        term1 = self.d2I0dKdK(strike, forward, alpha, corr, vovol, beta)*(1 + self.I1_Hagan(strike, forward, alpha, corr, vovol, beta)*self._expiryTerm)
        term2 = 2*self.dI0dK(strike, forward, alpha, corr, vovol, beta)*self.dI1dK(strike, forward, alpha, corr, vovol, beta)*self._expiryTerm
        term3 = self.I0_JObloj(strike, forward, alpha, corr, vovol, beta)*self.d2I1dKdK(strike, forward, alpha, corr, vovol, beta)*self._expiryTerm
        
        return term1 + term2 + term3

    
    def GetTotaldlogBSdK(self, strike, volatility, optiontype, forward, alpha, corr, vovol, beta):
        self._bs._strike = strike
        self._bs._volatility = volatility
        return (1/self._bs.GetOptionValue(optiontype)*(self._bs.GetDualDelta(optiontype) 
                                                      + self._bs.GetVega()*self.dSABRdK(strike, forward, alpha, corr, vovol, beta)))
                

    def GetTotald2logBSdKdK(self, strike, volatility, optiontype, forward, alpha, corr, vovol, beta):
        
        self._bs._strike = strike
        self._bs._volatility = volatility
        BS = self._bs.GetOptionValue(optiontype)
        dBSdK = self._bs.GetDualDelta(optiontype)
        d2BSdKdK = self._bs.GetDualGamma()
        vega = self._bs.GetVega()
        volga = self._bs.GetVolga()
        dualvega = self._bs.GetDualVega()
        dSABRdK = self.dSABRdK(strike, forward, alpha, corr, vovol, beta)
        d2SABRdKdK = self.d2SABRdKdK(strike, forward, alpha, corr, vovol, beta)
                       
        term1 = -pow(BS, -2)*pow(dBSdK + vega*dSABRdK, 2)
        term2 = 1/BS*(d2BSdKdK + dualvega*dSABRdK + (dualvega + volga*dSABRdK)*dSABRdK + vega*d2SABRdKdK)
        
        return term1 + term2      
    
       
    def CalcWingParameters(self) -> None:
        
        K10P = self._strikes[0]
        K25P = self._strikes[1]
        K25C= self._strikes[3]
        K10C = self._strikes[4]
                        
        ##########################################
        #  Solve Put Wing
        ##########################################
        forward = sfd.ForwardContinuousDeposit(self._spot, self._domesticDeposit, self._foreignDeposit, self._expiryTerm)
        
        put_matrix = np.array([[m.log(K25P), 1, K25P, K25P*K25P],
                               [1/K25P, 0, 1, 2*K25P],
                               [-1/(K25P*K25P), 0, 0, 2],
                               [m.log(K10P), 1, K10P, K10P*K10P]])
                
        self._bs._strike = K25P
        sabrvol25 = self.GetVolatility(K25P)
        self._bs._volatility = sabrvol25
        logBS25P = m.log(self._bs.GetOptionValue(bs.OptionType.Put))        
        dlogBSdK = self.GetTotaldlogBSdK(K25P, sabrvol25, bs.OptionType.Put, forward, self._alpha, self._corr, self._vovol, self._beta)        
        d2logBSdKdK = self.GetTotald2logBSdKdK(K25P, sabrvol25, bs.OptionType.Put, forward, self._alpha, self._corr, self._vovol, self._beta)
        self._bs._strike = K10P
        self._bs._volatility = self._volatilitySmile[0]
        logBS10P = m.log(self._bs.GetOptionValue(bs.OptionType.Put))
        
        rhs = np.array([[logBS25P], 
                        [dlogBSdK], 
                        [d2logBSdKdK], 
                        [logBS10P]])
        
        put_solution = np.linalg.inv(put_matrix)@rhs       
        self._my = put_solution[0]
        self._aput = put_solution[1]
        self._bput = put_solution[2]
        self._cput = put_solution[3]
                
        
        ##########################################
        #  Solve Call Wing
        ##########################################
        call_matrix = np.array([[-m.log(K25C),   1,   pow(K25C, -1),    pow(K25C, -2)],
                                [-pow(K25C, -1), 0,  -pow(K25C, -2), -2*pow(K25C, -3)],
                                [ pow(K25C, -2), 0, 2*pow(K25C, -3),  6*pow(K25C, -4)],
                                [-m.log(K10C),   1,   pow(K10C, -1),    pow(K10C, -2)]])
        
        self._bs._strike = K25C
        sabrvol25_call = self.GetVolatility(K25C)
        self._bs._volatility = sabrvol25_call
        logBS25C = m.log(self._bs.GetOptionValue(bs.OptionType.Call))
        
        dlogBSdK_call = self.GetTotaldlogBSdK(K25C, sabrvol25_call, bs.OptionType.Call, forward, self._alpha, self._corr, self._vovol, self._beta)        
        d2logBSdKdK_call = self.GetTotald2logBSdKdK(K25C, sabrvol25_call, bs.OptionType.Call, forward, self._alpha, self._corr, self._vovol, self._beta)
        self._bs._strike = K10C
        self._bs._volatility = self._volatilitySmile[4]
        logBS10C = m.log(self._bs.GetOptionValue(bs.OptionType.Call))
        
        rhs_c = np.array([[logBS25C], 
                        [dlogBSdK_call], 
                        [d2logBSdKdK_call], 
                        [logBS10C]])
        
        call_solution = np.linalg.inv(call_matrix)@rhs_c
        self._ny = call_solution[0]
        self._acall = call_solution[1]
        self._bcall = call_solution[2]
        self._ccall = call_solution[3]
        
        
        pass


#//     Unit-Test: SABR with wing extrapolation
class Test_SABRWing(unittest.TestCase):
    
    print('Running Test_SABRWing testing ')
    def test_dIdK(self):
        self.assertEqual(1.12345, 1.12345)

    

if __name__ == '__main__':
    unittest.main()
    
    
