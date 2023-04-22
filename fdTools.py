# Tools applied for running finite difference engines

import unittest

"""
//	-----------------------------------   THOMAS ALGORITHM   ----------------------------------- //
//		Given sub diagonal a[M-1], diagonal b[M] and super diagonal c[M-1] in tridiagonal 
//		matrix, Q[a,b,c].
//		Righthand side r[M]. Solution to Q*a = r  =>  a = Q^-1 * r
//
//      NB!: a, b and c need to be input as decimal numbers. If a element is 2 and not 2.0 
//           Python will interpret the number to be an integer. 1/2 = 0 instead of 1.0/2.0 = 0.5
// 
//      Numerical Partial Differential Equations, J.W. Thomas, Springer-Verlag
//	-----------------------------------   ----------------   ----------------------------------- //
"""
def Thomas(a, b, c, r, M):

    c[0] = c[0]/b[0]
    r[0] = r[0]/b[0]
    
    for i in range(1, M-1, 1):
    
        b[i] = b[i] - a[i-1]*c[i-1]
        r[i] = r[i] - a[i-1]*r[i-1]       
                
        c[i] = c[i]/b[i]
        r[i] = float(r[i]/b[i])        
           
    b[M-1] = b[M-1] - a[M-2]*c[M-2]
    r[M-1] = r[M-1] - a[M-2]*r[M-2]
    
    #  Finding the solution vector going back up
    r[M-1] = r[M-1]/b[M-1]
    
    for i in range(M-2, -1, -1):
        r[i] = r[i] - c[i]*r[i+1]    
        
    return r




#// Unit-Test: fdTools:
import numpy as np
class TestfdTools(unittest.TestCase):
    
    def test_ThomasAlgorithm(self):
        r = np.array([-1.0, -4.0, -7.0])
        a = np.array([2.0, 4.0])
        b = np.array([1.0, 4.0, 1.0])
        c = np.array([1.0, -1.0])
        r_ = Thomas(a, b, c, r, 3)
        self.assertEqual((r_ == np.array([0.5, -3.0/2.0, -1.0])).all(), True)

    
if __name__ == '__main__':
    unittest.main()
