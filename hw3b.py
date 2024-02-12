#region imports
from hw2a import Simpson2
import math
#endregion

def tPDF(m):
    """Returns the t distribution probability density function, without the Km constant"""
    f1 = lambda x: (1 + ((x ** 2) / m)) ** -((m + 1) / 2) # p1071
    return f1


def Km(m):
    """Returns the Km constant for the final equation"""
    km = math.gamma((m / 2) + 0.5) / (math.sqrt(m * math.pi) * math.gamma(m / 2))
    return km
def  mash(m, z, z0):
    """Returns the probability using simspons integration of tPDF and multiplies it by the Km constant"""
    mashup = Km(m) * Simpson2(tPDF(m), z,z0, 100)
    return mashup

def main():
    """Main function, asks for inputs and calculates the probability using the mash function. Prints results to the CLI"""
    getout = False
    while getout is False:
        m = input('Enter a Degree of Freedom: ') # asking input for the degrees of freedom
        z = float(input('Enter a z value: ')) # asking input for the z value
        z0 = z - (z*50)                    # calculates the lower limit of integration
        print("Probability: ",round(mash(float(m), z, z0), 3))
        getout = input('Do you want to leave Y/N: ') == ("Y" or "y")
    return print("Ended")

if __name__ == '__main__':
    main()




