"""
Examples for running the anchor localization module on known data 
"""
import os,sys,inspect
import scipy.io as sio
curdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
rootdir = os.path.dirname(curdir)
modulesdDir = os.path.join(rootdir, 'modules')
dataDir = os.path.join(rootdir, 'data')
exampleDataDir = os.path.join(curdir, 'example-data')
sys.path.insert(0, modulesdDir)

from lps_classes import Solver, Solution, Settings
import numpy as np

if __name__ == "__main__":

    # Load data
    testData = sio.loadmat(os.path.join(exampleDataDir,'example_no_outlier_data'))

    # Initialize the solution post RANSAC
    solver = Solver()
    initial_solution = Solution()
    initial_solution.Bhat = testData['Bhat']
    initial_solution.D = testData['D']
    initial_solution.rows = testData['rows'] - 1
    initial_solution.columns = testData['columns'] - 1
    initial_solution.d = testData['d']
    initial_solution.inliers = testData['inliers']
    
    settings = Settings()
    # Call solver, but only the optimization part
    r, s = solver.find_optimal_solution(initial_solution, settings)
    
    