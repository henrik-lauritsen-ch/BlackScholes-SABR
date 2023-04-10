
import unittest
from Utility import Test_Utility
from SABRWing import Test_SABRWing
from VolatilitySurface import Test_VolSurface
from StrikeFromDelta import Test_StrikeFromDelta
from BlackScholes import TestBSMethods


def load_tests(loader, tests, pattern):
# Load Test Classes

    suite = unittest.TestSuite()
        
    suite.addTests(loader.loadTestsFromModule(Test_StrikeFromDelta()))
    suite.addTests(loader.loadTestsFromModule(Test_Utility()))
    suite.addTests(loader.loadTestsFromModule(Test_VolSurface()))
    suite.addTests(loader.loadTestsFromModule(Test_SABRWing()))
    suite.addTests(loader.loadTestsFromModule(TestBSMethods()))    
    
    return suite

if __name__ == '__main__':
    unittest.main()
    
    
    