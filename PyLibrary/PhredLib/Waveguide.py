# For calculating things related to waveguides
import math

# Speed of light in vacuum
C = 2.997925e8

class Rect:
    """Parameters relating to a particular mode in a cerain sized
air filled rectangular waveguide """

    # Waveguide size
    a = 0.02286
    b = 0.01016

    # Frequency
    f = 10e9

    # Mode
    m = 1
    n = 0
    
    def omega(self):
        "Angular frequency"
        return 2 * math.pi * self.f

    def k(self):
        "Wave number"
        return self.omega() / C

    def kc(self):
        "Cutoff wavenumber"
        return math.sqrt( (self.m * math.pi / self.a) ** 2 + \
                          (self.n * math.pi / self.b) ** 2 )

    def beta(self):
        "Propagation constant"
        return math.sqrt( self.k() ** 2 - self.kc() ** 2 )
    
    def lambda_c(self):
        "Cut off frequency of wave guide"
        return 2 * math.pi / self.kc()

    def lambda_g(self):
        "Guided wavelength"
        return 2 * math.pi / self.beta()

    def vp(self):
        "Phase velocity"
        return self.omega() / self.beta()

    
    
