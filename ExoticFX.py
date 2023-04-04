# Exotic FX Options under the Garman-Kohlhagen model - First gegneration FX exotics

import BlackScholes as bs
import unittest

class Exotic(bs.GarmanKohlhagen):
    def __init__(self, spot, strike, expiryTerm, depositDomestic, depositForeign, volatility, barrier, lowerBarrier=0):
        super().__init__(spot, strike, expiryTerm, depositDomestic, depositForeign, volatility)

        self._barrier = barrier
        self.lowerBarrier = lowerBarrier







#//     Unit-Test: FX Exotics
class Test_Exotics(unittest.TestCase):
    
    def test_one(self):
        self.assertEqual(1.123456, 1.123456)
        

if __name__ == '__main__':
   unittest.main()