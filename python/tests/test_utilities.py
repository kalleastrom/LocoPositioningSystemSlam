import unittest
import numpy as np
import scipy.io as sio
from numpy.linalg import eig

import os,sys,inspect

curdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
rootdir = os.path.dirname(curdir)
modulesdDir = os.path.join(rootdir, 'modules')
dataDir = os.path.join(rootdir, 'data')
testDir = os.path.join(rootdir, 'tests')
sys.path.insert(0,modulesdDir)

import lps_utilities as util
from lps_classes import Settings

if sys.version_info < (3, 3):
    from mock import MagicMock
else:
    from unittest.mock import MagicMock

class TestUtilities(unittest.TestCase):

    def setUp(self):
        # Loads data from any data file
        self.data = 1

    def test_safe_cholesky_factorization(self):
        # Fixture
        N = np.random.randint(5,10)
        H = np.random.rand(N,N)
        H_POS_DEF = H + np.transpose(H)
        H_NEG_DEF = H_POS_DEF - max(np.real(eig(H_POS_DEF)[0])) * np.eye(N)
        
        # Test
        Hp, Lp = util.safe_cholesky_factorization(H_POS_DEF)
        Hn, Ln = util.safe_cholesky_factorization(H_NEG_DEF)
        
        # Assert
        self.assertAlmostEqual(Hp.all(), np.dot(np.transpose(Lp), Lp).all())
        self.assertAlmostEqual(Hn.all(), np.dot(np.transpose(Ln), Ln).all())

    def test_updatexy(self):
        # Fixture
        Nx = np.random.randint(5,6)
        Ny = np.random.randint(5,6)
        x = np.random.rand(3,Nx)
        y = np.random.rand(3,Ny)
        dz = np.random.rand(3*(Nx + Ny),1)
        dzx = dz[0:(3 * Nx), :]
        dzy = dz[(3 * Nx):, :]
        xnew = x + np.reshape(dzx, (3, Nx))
        ynew = y + np.reshape(dzy, (3, Ny))

        # Test
        xhat, yhat = util.updatexy(x, y, dz)

        # Assert
        self.assertAlmostEqual(xnew.all(), xhat.all())
        self.assertAlmostEqual(ynew.all(), yhat.all())
    
    def test_calcresandjac(self):
        # Fixture
        testData = sio.loadmat(os.path.join(testDir,'data_calcresandjac'))
        D = testData['D']
        I = testData['I'] - 1 # Different indexation in Python and Matlab
        J = testData['J'] - 1 # Different indexation in Python and Matlab
        x = testData['x']
        y = testData['y']

        # Assert
        res, jac = util.calcresandjac(D, I, J, x, y)

        # Test
        self.assertAlmostEqual(res.all(), testData['res'].all())
        self.assertAlmostEqual(jac.todense().all(), testData['jac'].all())

    def test_bundletoa(self):
        # Fixture
        testData = sio.loadmat(os.path.join(testDir,'data_bundletoa'))
        D = testData['D']
        I = testData['I'] - 1 # Different indexation in Python and Matlab
        J = testData['J'] - 1 # Different indexation in Python and Matlab
        x = testData['x']
        y = testData['y']
        
        settings = Settings()
        settings.bundle.numberOfIterations = 30
        settings.bundle.counterLimit = 50
        settings.bundle.numericalLimit = 1e-4

        # Assert
        xopt, yopt, res, jac = util.bundletoa(D, I, J, x, y, settings)

        # Test
        self.assertAlmostEqual(xopt.all(), testData['xopt'].all())
        self.assertAlmostEqual(yopt.all(), testData['yopt'].all())
        self.assertAlmostEqual(res.all(), testData['res'].all())
        self.assertAlmostEqual(jac.todense().all(), testData['jac'].todense().all())

if __name__ == '__main__':

    # Assert
    unittest.main()