# Some handy functions for plotting

import biggles
import Numeric

def plot_signal(signal, tsteps, dt):
    """Plot an input signal. Use this to see how your input signal looks.

signal is the SouceFunction object to plot
tsteps is the number of time steps
dt is the time step size
"""
    
    time = Numeric.arange(0, dt * tsteps, dt)
    y = [signal.signal_function(x) for x in time]
    
    p = biggles.FramedPlot()
    p.xlabel = "Time (sec)"
    p.ylabel = "Signal"
    p.add( biggles.Curve(time, y, color="red") )
    p.show()  
