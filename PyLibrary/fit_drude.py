# Try to fit the Drude model to the data from Johnson and Christy

from Numeric import *
from pylab import *
import pyx

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
    
def write_plot_gold():

    # Chang et. al. data:
    eps_inf = 11.4577
    wp = 9.4027 / H_ev * 2 * pi
    v = 0.08314 / H_ev

    # My fit using lqsnonlin in Matlab:
    #eps_inf = 8.09179
    #wp = 8.79336 / H_ev * 2 * pi
    #v = 0.452878 / H_ev

    jc_data = load_data('Johnson_Christy_gold.data')

    jc_n = jc_data[:,1] + jc_data[:,2] * 1j
    jc_eps = jc_n**2

    f = jc_data[:,0] / H_ev
    wl = C / f
    wld = wl * 1e9

    eps = drude(f * 2 * pi, eps_inf, wp, v)

    g = pyx.graph.graphxy(width=9.6, height=7.2, #width=12, height=9,
                          key=pyx.graph.key.key(), #hdist=2),
                          x=pyx.graph.axis.linear(min=300, max=1000,
                                                  title='Wavelength (nm)'),
                          y=pyx.graph.axis.linear(min=-50, max=10,
                                                  title='Real permittivity'),
                          y2=pyx.graph.axis.linear(min=-1, max=8,
                                                   title='Imaginary permittivity'))

    data_real1 = transpose(array((wld, jc_eps.real)))
    data_real2 = transpose(array((wld, eps.real)))
    data_imag1 = transpose(array((wld, jc_eps.imag)))
    data_imag2 = transpose(array((wld, eps.imag)))
    
    g.plot(pyx.graph.data.list(data_real1.tolist(), x=1, y=2,
                               title='J \& C $\Re\{\epsilon\}$'),
           [pyx.graph.style.line([pyx.style.linewidth.Thick,
                                  #pyx.style.linewidth.Thin,
                                  pyx.style.linestyle.dashed,
                                  pyx.color.rgb.blue])])
    g.plot(pyx.graph.data.list(data_real2.tolist(), x=1, y=2,
                               title='Drude $\Re\{\epsilon\}$'),
           [pyx.graph.style.line([pyx.style.linewidth.Thick,
                                  #pyx.style.linewidth.Thin,
                                  pyx.style.linestyle.solid,
                                  pyx.color.rgb.blue])])

    g.plot(pyx.graph.data.list(data_imag1.tolist(), x=1, y2=2,
                               title='J \& C $\Im\{\epsilon\}$'),
           [pyx.graph.style.line([pyx.style.linewidth.Thick,
                                  #pyx.style.linewidth.Thin,
                                  pyx.style.linestyle.dashed,
                                  pyx.color.rgb.green])])
    g.plot(pyx.graph.data.list(data_imag2.tolist(), x=1, y2=2,
                               title='Drude $\Im\{\epsilon\}$'),
           [pyx.graph.style.line([pyx.style.linewidth.Thin,
                                  pyx.style.linestyle.solid,
                                  pyx.color.rgb.green])])    
    
    g.writeEPSfile("gold_permittivity.eps") 

def write_plot_silver():
    # Hand calculations
    eps_inf = 4.15
    wp = 2.17e15*2*pi
    v = 32.258e12

    # My fit using lqsnonlin in Matlab:
    #eps_inf = 3.87717
    #wp = 9.17494 / H_ev * 2 * pi
    #v = 0.122595 / H_ev
    
    jc_data = load_data('Johnson_Christy_silver.data')

    jc_n = jc_data[:,1] + jc_data[:,2] * 1j
    jc_eps = jc_n**2

    f = jc_data[:,0] / H_ev
    wl = C / f
    wld = wl * 1e9

    eps = drude(f * 2 * pi, eps_inf, wp, v)

    g = pyx.graph.graphxy(width=9.6, height=7.2, #width=12, height=9,
                          key=pyx.graph.key.key(), #hdist=2),
                          x=pyx.graph.axis.linear(min=300, max=1000, #min=0, max=2000,
                                                  title='Wavelength (nm)'),
                          y=pyx.graph.axis.linear(min=-60, max=10,
                                                  title='Real permittivity'),
                          y2=pyx.graph.axis.linear(min=-1, max=3,
                                                   title='Imaginary permittivity'))

    data_real1 = transpose(array((wld, jc_eps.real)))
    data_real2 = transpose(array((wld, eps.real)))
    data_imag1 = transpose(array((wld, jc_eps.imag)))
    data_imag2 = transpose(array((wld, eps.imag)))
    
    g.plot(pyx.graph.data.list(data_real1.tolist(), x=1, y=2,
                               title='J \& C $\Re\{\epsilon\}$'),
           [pyx.graph.style.line([pyx.style.linewidth.Thick,
                                  #pyx.style.linewidth.Thin,
                                  pyx.style.linestyle.dashed,
                                  pyx.color.rgb.blue])])
    g.plot(pyx.graph.data.list(data_real2.tolist(), x=1, y=2,
                               title='Drude $\Re\{\epsilon\}$'),
           [pyx.graph.style.line([pyx.style.linewidth.Thick,
                                  #pyx.style.linewidth.Thin,
                                  pyx.style.linestyle.solid,
                                  pyx.color.rgb.blue])])

    g.plot(pyx.graph.data.list(data_imag1.tolist(), x=1, y2=2,
                               title='J \& C $\Im\{\epsilon\}$'),
           [pyx.graph.style.line([pyx.style.linewidth.Thick,
                                  #pyx.style.linewidth.Thin,
                                  pyx.style.linestyle.dashed,
                                  pyx.color.rgb.green])])
    g.plot(pyx.graph.data.list(data_imag2.tolist(), x=1, y2=2,
                               title='Drude $\Im\{\epsilon\}$'),
           [pyx.graph.style.line([pyx.style.linewidth.Thick,
                                  #pyx.style.linewidth.Thin,
                                  pyx.style.linestyle.solid,
                                  pyx.color.rgb.green])])    
    
    g.writeEPSfile("silver_permittivity.eps") 


def gold():
    # My fit using lqsnonlin in Matlab:
    eps_inf = 8.09179
    wp = 8.79336 / H_ev * 2 * pi
    vc = 0.452878 / H_ev

    jc_au = load_data('Johnson_Christy_gold.data')

    plot_drude(eps_inf, wp, v, jc_au)

def silver():
    # My fit using lqsnonlin in Matlab:
    eps_inf = 3.87717
    wp = 9.17494 / H_ev * 2 * pi
    vc = 0.122595 / H_ev
    
    print "Silver plasma freq: %g eV" % (wp / (2 * pi) * H_ev)
    print "Silver collision frequency: %g eV" % (v * H_ev)

    jc_ag = load_data('Johnson_Christy_silver.data')
    
    plot_drude(eps_inf, wp, v, jc_ag)
    
if (__name__ == "__main__"):
    
    
    write_plot_gold()
    write_plot_silver()
