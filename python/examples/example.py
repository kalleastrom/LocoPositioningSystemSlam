"""
Example when running the anchor localization module on known data 
"""
import os,sys,inspect

curdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
rootdir = os.path.dirname(curdir)
modulesdDir = os.path.join(rootdir, 'modules')
dataDir = os.path.join(rootdir, 'data')
sys.path.insert(0,modulesdDir)

from locolib import AnchorLocalizer
import numpy as np

if __name__ == "__main__":
    
    # Required inputs
    numberOfAnchors = 6
    filename = 'data4.json'
    print 10*'-'
    
    # Initialize solver
    solver = AnchorLocalizer(dataDir)
    
    # Prints the current state of the solver without loaded data
    print solver.verbose()
    print 10*'-'
    
    # Load data - checking that it is syntactically correct
    solver.load(filename)

    # Prints the current state of the solver with loaded data and no solution
    print solver.verbose()
    print 10*'-'
    
    # Call solver
    solver()
    
    # Print statistics
    print solver.verbose()
    print 2*'-'

    
    solver.initialPositions = np.random.rand(3,6)
    solver.finalPositions = np.random.rand(3,6)
    solver.finalPositionDeviations = 0.2*np.random.rand(3,6)
    solver.plot(['ranges','anchors'])