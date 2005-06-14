# Try to fit the Drude model to the data from Johnson and Christy

from Numeric import *
from pylab import *

H_ev = 4.13566743e-15 # eV s
H_j = 6.6260693e-34 # J s
C = 2.979e8 # m/s

def load_data(filename):
    f = open(filename, 'r')

    data = f.read().split('\n')
    f.close()

    arr = zeros((len(data), 3), typecode='f')

    j = 0
    for line in data:
        try:
            items = [float(x) for x in line.split(" ")]
        except ValueError:
            continue
        
        if (len(items) != 3):
            continue
        
        for i in range(3):
            arr[j, i] = items[i]

        j = j + 1

    return arr[0:j,:]

def drude(w, eps_inf, wp, v):
    return eps_inf - (wp**2) / ((w**2) + 1j * w * v)


def plot_drude(eps_inf, wp, v, jc_data):
    jc_n = jc_data[:,1] + jc_data[:,2] * 1j
    jc_eps = jc_n**2

    f = jc_data[:,0] / H_ev
    wl = C / f
    wld = wl * 1e9

    eps = drude(f * 2 * pi, eps_inf, wp, v)

    plot(wld, jc_eps.real, wld, eps.real)
    legend(('JC real', 'Test real'))

    figure()
    plot(wld, jc_eps.imag, wld, eps.imag)
    legend(('JC Imag', 'Test imag'))
    show()
    

def gold():
    eps_inf = 11.4577
    wp = 9.4027 / H_ev * 2 * pi
    v = 0.08314 / H_ev

    jc_au = load_data('Johnson_Christy_gold.data')

    plot_drude(eps_inf, wp, v, jc_au)

def silver():
    eps_inf = 1
    wp = 2.17e15*2*pi
    v = 32.258e12

    jc_ag = load_data('Johnson_Christy_silver.data')

    plot_drude(eps_inf, wp, v, jc_ag)
    
    
if (__name__ == "__main__"):
    gold()
