"""
Examples for running the dicrete time dynamical models open loop 
"""
import os,sys,inspect
from math import ceil, sin
import matplotlib.pylab as plt

curdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
rootdir = os.path.dirname(curdir)
modulesdDir = os.path.join(rootdir, 'modules')
sys.path.insert(0,modulesdDir)

from lps_classes import TrippleIntegrator

import numpy as np

if __name__ == "__main__":

    # Initialize a 3D tripple integrator, discretized at 0.1 s
    timeSimulation = 20.0
    timeStep = 0.1
    TI = TrippleIntegrator(dt=timeStep)
    
    t, u, x, y = [], [], [], []
    
    # Control signal sequence
    controlSignal = lambda t: np.array([[1.0*sin(t)], [-0.3*sin(3.0*t)], [0.5*sin(t + 1.0)]])
    
    # Simulate over 10 seconds
    for ii in range(int(ceil(timeSimulation/timeStep))):
        time = ii * timeStep
        csig = controlSignal(time)
        t = t + [time]
        u = u + np.transpose(csig).tolist()
        state, measurement = TI(csig, timeStep)
        x = x + np.transpose(state).tolist()
        y = y + measurement.tolist()

    # Plot results
    plt.figure(1)
    plt.step(t, u)
    plt.title('Control signals')
    plt.xlabel('Time [$s$]')
    plt.ylabel('Input signal [$m/s^2$]')

    plt.figure(2)
    plt.step(t, u)
    plt.title('Measurements')
    plt.xlabel('Time [$s$]')
    plt.ylabel('Measurement [$m/s$]')

    plt.figure(3)
    x = np.array(x)
    plt.subplot(311)
    plt.step(t, x[:,0:3])
    plt.title('Position')
    plt.xlabel('Time [$s$]')
    plt.ylabel('Measurement [$m$]')
    plt.subplot(312)
    plt.step(t, x[:,3:6])
    plt.title('Velocity')
    plt.xlabel('Time [$s$]')
    plt.ylabel('Measurement [$m/s$]')
    plt.subplot(313)
    plt.step(t, x[:,6:9])
    plt.title('Acceleration')
    plt.xlabel('Time [$s$]')
    plt.ylabel('Measurement [$m/s^2$]')
    plt.tight_layout(pad=0.01)