"""
 FX Vanilla option tools: normal distribution, interpolation, solvers and other

"""
import math
import numpy as np
import unittest
import sys

class Constants:
    
    def __init__(self) -> None:
        pass
    
    _Pi = 3.1415926535897932
    _SqrtTwo = 1.4142135623730951
    _OneOverRootTwoPi = 0.398942280401433



# This method maps a real numer y to the interval [a, b] using the inverse of 
# y = (x-(a+b)/2)/((x-a)(x-b)). y is a bijektiv mapping of [a, b] into R.
# The property of y:
# y((a+b)/2) = 0.0 
# y ->  infinity for x -> a+ 
# y -> -infinity for x -> b-
# The inverse of y is given as
# x = (y*(a+b+1/y) + sqrt(y^2*(a+b+1/y)^2 - 4*y*(a*b*y+(a+b)/2)))/2*y
def RealAxisToIntervalAB(y, a, b):

    a1 = y
    b1 = -(y * a + y * b + 1)
    c1 = a * b * y + (a + b) / 2

    if (y == 0.0):
        retval = (a + b) / 2
    else:
        retval = SolutionToSecondDegreePolynomial(a1, b1, c1, False)

    return retval


# This method maps a real numer from the interval [a, b] into R using 
# y = (x-(a+b)/2)/((x-a)(x-b)). y is bijective
def IntervalABToRealAxis(x, a, b):

    if (x == a):
        retval = sys.float_info.max
    elif (x == b):
        retval = sys.float_info.min
    else:
        retval = (x - (a + b) / 2) / ((x - a) * (x - b))

    return retval

# This method returne the solution to the second order polynomial equation:
# a*x^2 + b*x + c = 0.
# Solution: x = (-b +/- sqrt(b^2 - 4*a*c))/(2*a)
# If plussoltion is TRUE then solution x = (-b + sqrt(b^2 - 4*a*c))/(2*a) is returned
def SolutionToSecondDegreePolynomial(a, b, c, plussolution = True):

    if ((b * b - 4 * a * c) >= 0.0):
        sqrtbb4ac = math.sqrt(b * b - 4 * a * c)
    else:
        raise ValueError('No solution to second degree polynomial')

    if (plussolution == True):
        retval = (-b + sqrtbb4ac) / (2 * a)
    else:
        retval = (-b - sqrtbb4ac) / (2 * a)

    return retval


def bisection(func, x1, x2, x_accuracy, JMAX = 40):
    # Method taken from Numerical Recipes in C
    # https://web.stanford.edu/~doubleh/eco270/bisection
    f = func(x1)
    fmid = func(x2)
        
    if (f*fmid>0):
        raise ValueError('No root in interval from min_x to max_x')
    
    rtb = 0
    if f<0:
        dx, rtb = x2 - x1, x1
    else:
        dx, rtb = x1 - x2, x2        
           
    i = 1   
    while (i <= JMAX):
        dx *= 0.5
        xmid = rtb + dx
        fmid = func(xmid)
                
        if (fmid <= 0.0):
            rtb = xmid             
    
        if (abs(dx) < x_accuracy or fmid == 0.0): 
            break
        
        i += 1
        
    return rtb

    
class norm:    
    
    def __init__(self) -> None:
        pass

    
    def pdf(self, x):
        return Constants._OneOverRootTwoPi * math.exp(-x * x / 2)


    def InfErf(self, x):
            
        factor = 1.12837916709551
        sign, sgnflag = 1, 0

        if (x < 0):
            sign = -1
            x = -x
            sgnflag = 1

        if (x == 0):
            ret = 0
        elif (x < 3.6):	 #apply series solution
            s = 1
            sum = 1
            xN = 0
            x2Nm1 = -1
            x2Np1 = 1
            x2 = x * x

            if (x < 0.5):
                n = 9
            else:
                tmp = 46 * (x - 0.5) / 3.1
                if (tmp > 0):
                    n = int(math.ceil(tmp)) + 9
                else:
                    n = int(math.ceil(tmp)) + 9

            # int N=(x<0.5) ? 9:((int)(ceil(46*(x-0.5)/3.1)+9))
            i = 1
            while (i <= n):

                xN = xN + 1
                x2Nm1 = x2Nm1 + 2
                x2Np1 = x2Np1 + 2
                s = s * (-x2 * x2Nm1 / (xN * x2Np1))
                sum = sum + s
                i = i + 1

            ret = sum * x * factor
                
        elif (x <= 5.4):  # apply asymptotic solution           
            s = 1
            sum = 1
            xN = 0
            x2N = 0
            x2Nm1 = -1
            x2 = 4 * x * x

            i = 1
            while (i <= 14):        
                xN = xN + 1
                x2N = x2N + 2
                x2Nm1 = x2Nm1 + 2
                s = s * (-x2Nm1 * x2N / (x2 * xN))
                sum = sum + s
                i = i + 1
            
            ret = 1 - math.exp(-x * x) * factor * sum / (2 * x)
                
        else:
            ret = 1.0


        if (sgnflag == 1):
            ret = ret * sign

        return ret
        

    def cdf(self, x):       
        return self.cdfI(x)


    def cdfM(self, x):
        retval = 0    
        y = 1.0 / (1 + 0.2316419 * abs(x))
        z = 1.330274429 * math.pow(y, 5) - 1.821255978 * math.pow(y, 4) + 1.781477937 * math.pow(y, 3) - 0.356563782 * math.pow(y, 2) + 0.31938153 * y
        nu = 1.0 - z / math.sqrt(2 * Constants._Pi) * math.exp(-math.pow(x, 2) / 2)

        if (x > 0.0):
            retval = nu
        else:
            retval = 1 - nu

        return retval


    def cdfI(self, x):
        SqrtTwo = Constants._SqrtTwo
        tmp = 0.5 * (self.InfErf(x / SqrtTwo) + 1)
        return tmp


    def InverseCdf(self, x):
        """   Inverse normal algorithm developed by
            P. J. Acklam.  It is accurate to 1.15E-9.  Method checked
            against the Royal Statistical Society algorithm 241
            for single precision (Ppnd7).  The agreement for p = 0.01
            to p = 0.99 is no worse than 1.0E-7, the accuracy of Ppnd7.
        """
        a1 = -39.6968302866538
        a2 = 220.946098424521
        a3 = -275.928510446969
        a4 = 138.357751867269
        a5 = -30.6647980661472
        a6 = 2.50662827745924

        b1 = -54.4760987982241
        b2 = 161.585836858041
        b3 = -155.698979859887
        b4 = 66.8013118877197
        b5 = -13.2806815528857

        c1 = -7.78489400243029E-03
        c2 = -0.322396458041136
        c3 = -2.40075827716184
        c4 = -2.54973253934373
        c5 = 4.37466414146497
        c6 = 2.93816398269878

        d1 = 7.78469570904146E-03
        d2 = 0.32246712907004
        d3 = 2.445134137143
        d4 = 3.75440866190742

        p_low = 0.02425
        p_high = 1 - p_low

        retval = 0
        if (x < 0 or x > 1):
            return 1.0
        elif (x < p_low):            
            q = math.sqrt(-2 * math.log(x))
            retval = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /((((d1 * q + d2) * q + d3) * q + d4) * q + 1)
        elif (x <= p_high):
            q = x - 0.5
            r = q * q
            retval = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /(((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1)
        else:
            q = math.sqrt(-2 * math.log(1 - x))
            retval = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /((((d1 * q + d2) * q + d3) * q + d4) * q + 1)

        return retval


    def Moro(self, u):
        
        # Moro Inverse Cumulative Normal
        a = np.array ([2.50662823884, 
                    -18.61500062529, 
                    41.39119773534, 
                    -25.44106049637],
                    dtype=float
                    )

        b = np.array([-8.47351093090,
                    23.08336743743,
                    -21.06224101826,
                    3.13082909833],
                    dtype=float
                    )
        
        c = np.array([0.3374754822726147,
                    0.9761690190917186,
                    0.1607979714918209,
                    0.0276438810333863,
                    0.0038405729373609,
                    0.0003951896511919,
                    0.0000321767881768,
                    0.0000002888167364,
                    0.0000003960315187],
                    dtype=float
                    )

        x = u - 0.5
        r = 0
        
        if (abs(x) < 0.42):  # Beasley-Springer
            y = x * x
            r = x * (((a[3] * y + a[2]) * y + a[1]) * y + a[0]) /((((b[3] * y + b[2]) * y + b[1]) * y + b[0]) * y + 1.0)
                
        else: # Moro
            r = u
            if (x > 0.0):
                r = 1.0 - u

            r = math.log(-math.log(r))
            r = c[0] + r * (c[1] + r * (c[2] + r * (c[3] + r * (c[4] + r * (c[5] + r * (c[6] + r * (c[7] + r * c[8])))))))

            if (x < 0.0):
                r = -r

        return r


def FindIndex(x, data):
    """     Given Container data[num_cols], and given a value x, a value index
    is returned such that x is between data[index] and data[index+1]. data must
    be monotonic, either increasing or decreasing. index = -1 or index = num_cols
    is returned to indicate that x is out of range
    """
    num_cols = len(data)
    klo = 0
    khi = num_cols-1

    if (x < data[0]):
        index = -1
    elif (x > data[num_cols-1]):
        index = num_cols
    else:
        while(khi - klo > 1):

            k = int((khi+klo)/2)
            if (data[k] > x):
                khi = k
            else:
                klo = k
		
        index = klo
	
    return index


class Interpolation:
    
    def __init__(self, flatExtrapolation=False):
        self._flatExtrapolation = flatExtrapolation
    
    
    def GetInterpolatedValue(self, x, xs, ys):
        pass 
    
class PiecewiseLinearInterpolation(Interpolation):

    def __init__(self, flatExtrapolation=False) -> None:
        super().__init__(flatExtrapolation)    
    
    
    def LinearInterpolation(self, x, x1, x2, y1, y2):

        if (x2 - x1 != 0):
            y = y1 + (y2 - y1)/(x2 - x1)*(x - x1)
        else:
            y = y1

        return y


    def GetInterpolatedValue(self, x, xs, ys):
        
        nx = len(xs)
        
        if (x<xs[0]):
            if (self._flatExtrapolation==True):
                y = ys[0]
            else:
                y = self.LinearInterpolation(x, xs[0], xs[1], ys[0], ys[1])
        elif (x>xs[nx - 1]):
            if (self._flatExtrapolation==True):
                y = ys[nx - 1]
            else:
                y = self.LinearInterpolation(x, xs[nx - 2], xs[nx - 1], ys[nx - 2], ys[nx - 1])
        else:
             i = FindIndex(x, xs)                
             y = self.LinearInterpolation(x, xs[i], xs[i+1], ys[i], ys[i+1])        
        
        return y


class CubicSplineInterpolation(Interpolation):
    
    def __init__(self, flatExtrapolation=False):
        super().__init__(flatExtrapolation)
    
    def GetInterpolatedValue(self, x, xs, ys):
        # Method taken from Numerical Recipes in C
        # http://phys.uri.edu/nigh/NumRec/bookfpdf/f3-3.pdf
        
        nx = len(xs)        
        if (x < xs[0] and self._flatExtrapolation==True):
            y = ys[0]
        elif (x > xs[nx - 1] and self._flatExtrapolation==True):
            y = ys[nx - 1]
        else:
            y2s = self.spline(xs, ys)
            y = self.splint(x, xs, ys, y2s)
        
        return y
        

    def spline(self, xs, ys, yp1=0.99e31, ypn=0.99e31):

        numberOfX = len(xs)             
        u = np.zeros(numberOfX - 1)
        y2 = np.zeros(numberOfX)
        
        if (yp1 > 0.99e30):
            y2[0], u[0] = 0.0, 0.0
        else:
            y2[0] = -0.5
            u[0] = (3.0 / (xs[1] - xs[0])) * ((ys[1] - ys[0]) / (xs[1] - xs[0]) - yp1)
        
        for i in range(1, numberOfX-1, 1):
        
            sig = (xs[i] - xs[i - 1]) / (xs[i + 1] - xs[i - 1])
            p = sig * y2[i - 1] + 2.0
            y2[i] = (sig - 1.0) / p
            u[i] = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]) - (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1])
            u[i] = (6.0 * u[i] / (xs[i + 1] - xs[i - 1]) - sig * u[i - 1]) / p
        
        
        if (ypn > 0.99e30):
            qn, un = 0.0, 0.0
        else:
            qn = 0.5
            un = (3.0 / (xs[numberOfX - 1] - xs[numberOfX - 2])) * (ypn - (ys[numberOfX - 1] - ys[numberOfX - 2]) / (xs[numberOfX - 1] - xs[numberOfX - 2]))
        
        y2[numberOfX - 1] = (un - qn * u[numberOfX - 2]) / (qn * y2[numberOfX - 2] + 1.0)
        
        for k in range(numberOfX - 2, -1, -1):
            y2[k] = y2[k] * y2[k + 1] + u[k]
        
        return y2    


    def splint(self, x, xs, ys, y2s):

        y = 0
        numberofx = len(xs)

        if (x <= xs[0]):
            
            alpha = (ys[1] - ys[0]) / (xs[1] - xs[0]) - 1.0 / 6.0 * (xs[1] - xs[0]) * (2.0 * y2s[0] + y2s[1])
            y = ys[0] + alpha * (x - xs[0])
        
        elif (x >= xs[numberofx - 1]):
        
            alpha = ((ys[numberofx - 1] - ys[numberofx - 2]) / (xs[numberofx - 1] - xs[numberofx - 2]) +
                            1.0 / 6.0 * (xs[numberofx - 1] - xs[numberofx - 2]) * (y2s[numberofx - 2] + 2 * y2s[numberofx - 1]))
            y = ys[numberofx - 1] + alpha * (x - xs[numberofx - 1])
        
        else:    
            klo = 0
            khi = numberofx - 1
            while (khi - klo > 1):
            
                k = int((khi + klo) / 2)
                if (xs[k] > x):
                    khi = k
                else:
                    klo = k
            

            h = xs[khi] - xs[klo]
            if (h == 0.0):        
                raise ValueError("xs[i]!=xs[j] for all i,j where i!=j")
            
            a = (xs[khi] - x) / h
            b = (x - xs[klo]) / h
            y = a * ys[klo] + b * ys[khi] + ((a * a * a - a) * y2s[klo] + (b * b * b - b) * y2s[khi]) * (h * h) / 6.0
        
        return y



#//    Unit-Test: Functions of Utility
def g(x):
    return 2*x - 8.2468642

class Test_Utility(unittest.TestCase):
    
    def test_Bisection(self):
        self.assertEqual(round(bisection(g, 0, 7, 0.00000001) , 7), round(4.1234321, 7))

    def test_ICdf(self):
        self.assertEqual(round(norm().cdfI(2.134), 8), round(0.983578609590808, 8))
       
    def test_MCdf(self):
        self.assertEqual(round(norm().cdfM(3.0), 8), round(0.998650032777765, 8))

    def test_InverseCdf(self):
        self.assertEqual(round(norm().InverseCdf(0.9752), 8), round(1.96339753624734, 8))
    
    def test_Moro(self):
        self.assertEqual(round(norm().Moro(0.9752), 8), round(1.96339753624734, 8))        

    def test_CubicSplineInterpolation(self):
        xs = np.array([9.796265875871027, 10.067098505250692, 10.356101824110898, 10.687697656702378, 11.069582777590423])
        ys = np.array([0.09852,	0.09542, 0.0973, 0.10582, 0.11732])
            # Interpolation between 2 points
        cs = CubicSplineInterpolation(False)    
        
        self.assertEqual(round(cs.GetInterpolatedValue(10.7, xs, ys), 12), round(0.106188572556555, 12))
            # Linear extrapolation above right end point
        self.assertEqual(round(cs.GetInterpolatedValue(11.9, xs, ys), 12), round(0.14239424307285, 12))

    def test_PieceviseLinear(self):
        xa = np.array([2, 4.1, 7.3331, 9.998])
        ya = np.array([10, 7, 5.6, 11.1])
        pl = PiecewiseLinearInterpolation(True)
        self.assertEqual(pl.GetInterpolatedValue(1.1, xa, ya), 10)
        self.assertEqual(round(pl.GetInterpolatedValue(3.1, xa, ya), 14), 8.42857142857143)
        self.assertEqual(round(pl.GetInterpolatedValue(7.2, xa, ya), 14), 5.65763508706814)
        
        pl._flatExtrapolation = False
        self.assertEqual(round(pl.GetInterpolatedValue(21.1, xa, ya), 13), 34.0130548988705)

if __name__ == '__main__':
    unittest.main()
    