# region imports
from math import sqrt,exp,pi
# endregion

# region methods
def Probability(PDF, args, c, GT=True):
    '''
    Use Simpson rule to numerically integrate PDF.  Result is probability P(x<c).
    If I want probability x>c (GT=True), calculate 1-P
    :param PDF: the probability density function (a callback function)
    :param args: contains (mean, standard deviation) in that order
    :param c: value for calculating P(x<c)
    :param GT: boolean for deciding if we want P(x>c) (GT=True) or P(x<c) (GT=False)
    :return: P(x<c) or P(x>c)
    '''
    m, s = args  # unpack mean and standard deviation
    a = m-5*s  # compute left limit for integral as 5 standard deviations below (to left of) mean
    P = Simpson(PDF,args, a, c, npoints=100)  # use Simpson function to integrate PDF with 100 points

    return P if not GT else (1-P)  # return desired probability based on GT

def GNPDF(args):
    '''
    This is the Gaussian Normal Probability Density Function with parameters mean (mu), standard deviation (sigma).
    :param args: contains x, mean, standard deviation
    :return: f(x)
    '''
    x, mu, sigma = args  # unpack from args
    f = (1/(sigma*sqrt(2*pi)))*exp(-0.5*((x-mu)/sigma)**2)  # calculate f(x)  PEDMAS
    return f  # return f(x)

def Simpson(fcn, args, a, b, npoints = 21):
    '''
    This uses the simpson 1/3 rule for numerical integration, which uses quadratic Lagrange polynomials for interpolation.
    See Kryszig pp831-32 for derivation.
    The panel size is:  dX=abs(a-b)/(2*npoints)
    The integral is found by:  I=dX/3*(fcn(0)+fcn(npoints)+4*sum(fcn(j), 1 to n-1, odd)+2*(fcn(j), 2 to n-2, even))
    :param fcn: the callback function.  the function of x to integrate
    :param args: contains parameters for the Gaussian Normal distrubution (mean, standard deviation)
    :param a: left limit of integration (check a<b)
    :param b: right limit of integration (check a<b)
    :param npoints: number of points in integration
    :return: the value of the integral of fcn(x) between a,b
    '''
    mu, sig=args #unpack mu and sig
    area = 0  # initial value for the integral
    m = npoints  #staying consistent with book nomenclature
    n = 2*m  # ensure an even number of panels
    # if a,b are passed wrong, this checks to put them in correct order
    xL = min(a, b)  # lower limit of integration
    xR = max(a, b)  # upper limit of integration

    if xL == xR:
        return 0  # nothing to do, so return 0

    h = (xR-xL)/(2*m)  # calculate panel width
    x = xL
    area = (fcn((xL, mu,sig))+fcn((xR, mu,sig)))  # first calculate f(x) at xL and xR

    for j in range(1,2*m): # counts from 1 to 2*m-1
        x=j*h+xL # update the x position for evaluating the function
        if not j % 2 == 0: #the odds
            area += 4*fcn((x, mu,sig))
        else:  # the evens.
            area += 2*fcn((x, mu,sig))
    return (h/3.0)*area  # finally, return the value for the integral


def Simpson2(fcn, a, b, npoints = 21):
    '''
    This uses the simpson 1/3 rule for numerical integration, which uses quadratic Lagrange polynomials for interpolation.
    See Kryszig pp831-32 for derivation.
    The panel size is:  dX=abs(a-b)/(2*npoints)
    The integral is found by:  I=dX/3*(fcn(0)+fcn(npoints)+4*sum(fcn(j), 1 to n-1, odd)+2*(fcn(j), 2 to n-2, even))
    :param fcn: the callback function.  the function of x to integrate
    :param args: contains parameters for the Gaussian Normal distrubution (mean, standard deviation)
    :param a: left limit of integration (check a<b)
    :param b: right limit of integration (check a<b)
    :param npoints: number of points in integration
    :return: the value of the integral of fcn(x) between a,b
    '''
    area = 0  # initial value for the integral
    m = npoints  #staying consistent with book nomenclature
    n = 2*m  # ensure an even number of panels
    # if a,b are passed wrong, this checks to put them in correct order
    xL = min(a, b)  # lower limit of integration
    xR = max(a, b)  # upper limit of integration

    if xL == xR:
        return 0  # nothing to do, so return 0

    h = (xR-xL)/(2*m)  # calculate panel width
    x = xL
    area = (fcn((xL))+fcn((xR)))  # first calculate f(x) at xL and xR

    for j in range(1,2*m): # counts from 1 to 2*m-1
        x=j*h+xL # update the x position for evaluating the function
        if not j % 2 == 0: #the odds
            area += 4*fcn((x))
        else:  # the evens.
            area += 2*fcn((x))
    return (h/3.0)*area  # finally, return the value for the integral

def main():
    '''
    This main function calculates probabilities using the probability function I defined.
    Step 1:  define mean and standard deviation.
    Step 2: set c=1 for P(x<1)|N(0,1)) i.e., c=mu1 + 1*stdev
    Step 3: calculate probability by calling Probability with GNPDF as callback, (mu1, stdev1) as tuple argument, c1, and GT=False
    Step 4: repeat steps 1-3 for P(x>mu2+2*stdev2|N(175,3))
    Step 4: calculate P(-1<x<1|N(0,1)
    Step 5: output results to the console.
    :return: no return value from this function.
    '''

    # P(x<1|N(0,1))
    # step 1
    mu1=0
    stdev1 = 1

    # step 2
    c1 = 1

    # step 3
    #GNPDF is the Gaussian-Normal Probability Density Function (see def GNPDF)
    P1 = Probability(GNPDF,(mu1,stdev1),c1,GT=False)  # P(x<1.00|N(0,1))

    # P(x>181.0|N(175,3))
    # step 1 (part 2)
    mu2 = 175
    stdev2 = 3

    # step 2 (part 2)
    c2 = mu2+2*stdev2

    # step 3 (part 3)
    P2 = Probability(GNPDF, (mu2,stdev2),c2)

    # P(-2.00<x<2.00|N(0,1))
    # step 1 (part 3)
    mu3 = 0.0
    stdev3 = 1.0

    # step 2 (part 3)
    c3 = 2.0

    #step 3 (part 3)
    P3 = 1.0-2 * Probability(GNPDF, (mu3, stdev3), c3)

    # step 5
    print('P(x<{:0.2f}|N({:0.2f},{:0.2f}))={:0.3f}'.format(c1, mu1,stdev1,P1))
    print('P(x>{:0.2f}|N({:0.2f},{:0.2f}))={:0.3f}'.format(c2, mu2,stdev2,P2))
    print('P({:0.2f}<x<{:0.2f}|N({:0.2f},{:0.2f}))={:0.3f}'.format(-1.0*c3, c3, mu3, stdev3, P3))

# endregion

# region function call(s)
if __name__ == "__main__":
    main()
# endregion