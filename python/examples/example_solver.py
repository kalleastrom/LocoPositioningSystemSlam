"""
Examples for running the anchor localization module on known data 
"""
import os,sys,inspect

curdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
rootdir = os.path.dirname(curdir)
modulesdDir = os.path.join(rootdir, 'modules')
dataDir = os.path.join(rootdir, 'data')
sys.path.insert(0,modulesdDir)

from lps_classes import Solver
import numpy as np

if __name__ == "__main__":
    
    # Required inputs
    numberOfAnchors = 6
    filename = 'data4.mat'
    
    # Initialize solver
    b = 1
    solver = Solver(directory = dataDir)
    
    # Prints the current state of the solver without loaded data
    print solver.solution.verbose()
    
    # Load data - checking that it is syntactically correct
    solver.load(filename)

    # Prints the current state of the solver with loaded data and no solution
    print solver.solution.verbose()
    
    # Call solver
    #solver()
    
    # Print statistics
    print solver.solution.verbose()

    # These fields will be set in the call() function once unit tests pass
    solver.solution.initialPositions = np.random.rand(3,6)
    solver.solution.finalPositions = np.random.rand(3,6)
    solver.solution.finalPositionDeviations = 0.2*np.random.rand(3,6)
    solver.solution.plot(['anchors'])