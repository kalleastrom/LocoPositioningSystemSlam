"""
Examples for running the anchor localization module on known data 
"""
import os,sys,inspect
import scipy.io as sio
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D

curdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
rootdir = os.path.dirname(curdir)
modulesdDir = os.path.join(rootdir, 'modules')
dataDir = os.path.join(rootdir, 'data')
exampleDataDir = os.path.join(curdir, 'example-data')
sys.path.insert(0, modulesdDir)

from lps_classes import Solver, Solution, Settings

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
    
    # Create an instance of default settings
    settings = Settings()

    # Call solver, but only the optimization part
    r, s = solver.find_optimal_solution(initial_solution, settings)

    # Cisualize solution
    fig = plt.figure(1)
    ax = fig.gca(projection='3d')
    s = np.array([s[:,ii] for ii in range(s.shape[1]) if 2.0>max(s[:,ii])]).T
    ax.plot(s[0,:], s[1,:], s[2,:], color='r',alpha=0.5)
    ax.plot(r[0,:], r[1,:], r[2,:],'^k',markersize=10)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.legend(['Sender, s_i', 'Receiver, r_i'],loc=4)
    ax.set_title('Computed trajectory (S) and anchor positions (R) \n bundle toa_3D_bundle without smoothing')
    plt.show()
    