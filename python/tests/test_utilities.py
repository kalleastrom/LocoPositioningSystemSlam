import unittest
import numpy as np
import scipy.io as sio
from numpy.linalg import eig, inv
import matplotlib.pylab as plt

import os,sys,inspect
from scipy.sparse.linalg import lsqr
import scipy
import numpy

curdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
rootdir = os.path.dirname(curdir)
modulesdDir = os.path.join(rootdir, 'modules')
dataDir = os.path.join(rootdir, 'data')
testDir = os.path.join(rootdir, 'tests')
testDataDir = os.path.join(testDir, 'test-data')
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
        H = H + H.T
        H_POS_DEF = H + max(np.real(eig(H)[0])) * np.eye(N)
        H_NEG_DEF = H - max(np.real(eig(H)[0])) * np.eye(N)
        
        # Test
        Hp, Lp = util.safe_cholesky_factorization(H_POS_DEF)
        Hn, Ln = util.safe_cholesky_factorization(H_NEG_DEF)
        
        # Assert
        status = False
        try:
            # Positive definite matrix
            np.testing.assert_array_almost_equal(inv(Hp), Lp.dot(Lp.T), decimal=5)
            # Negative definite matrix
            np.testing.assert_array_almost_equal(inv(Hn), Ln.dot(Ln.T), decimal=5)
            # If matrix is positive definite, H should equal to Hp
            np.testing.assert_array_almost_equal(H_POS_DEF, Hp, decimal=5)
            status = True
        except AssertionError as e:
            print e.message
            pass
        self.assertTrue(status)

    def test_updatexy(self):
        # Fixture
        testData = sio.loadmat(os.path.join(testDataDir,'data_updatexy'))
        x = testData['x']
        y = testData['y']
        dz = testData['dz']

        # Test
        xhat, yhat = util.updatexy(x, y, dz)

        # Assert
        status = False
        try:
            np.testing.assert_array_almost_equal(xhat, testData['xny'], decimal=5)
            np.testing.assert_array_almost_equal(yhat, testData['yny'], decimal=5)
            status = True
        except AssertionError as e:
            print e.message
            pass
        self.assertTrue(status)
    
    def test_toa_normalize(self):
        # Fixture
        testData = sio.loadmat(os.path.join(testDataDir,'data_toa_normalize'))
        rin = testData['rin']
        sin = testData['sin']

        # Test        
        rout, sout = util.toa_normalise(rin, sin)

        # Assert
        status = False
        try:
            # Receiver positions
            np.testing.assert_array_almost_equal(rout, testData['rout'], decimal=5)
            # Sender positions
            np.testing.assert_array_almost_equal(sout, testData['sout'], decimal=5)
            status = True
        except AssertionError as e:
            print e.message
            pass
        self.assertTrue(status)

    def test_calcresandjac(self):
        # Fixture
        testData = sio.loadmat(os.path.join(testDataDir,'data_calcresandjac'))
        D = testData['D']
        I = testData['I'] - 1 # Different indexation in Python and Matlab
        J = testData['J'] - 1 # Different indexation in Python and Matlab
        x = testData['x']
        y = testData['y']

        # Test
        res, jac = util.calcresandjac(D, I, J, x, y)

        # Assert            
        status = False
        try:
            # Test the accuracy of the residual
            np.testing.assert_array_almost_equal(res, testData['res'], decimal=5)
            # Test the computed jacobian in full array form
            np.testing.assert_array_almost_equal(jac.toarray(), testData['jac'], decimal=5)
            status = True
        except AssertionError as e:
            print e.message
            pass
        self.assertTrue(status)

    def test_get_linear_constraints(self):
        # Fixture
        testData = sio.loadmat(os.path.join(testDataDir,'data_get_linear_constraints'))
        Bhat = testData['Bhat']
        xt = testData['xt']

        # Test
        H, b, cfm_linear = util.get_linear_constraints(Bhat, xt)

        print H.shape
        print b.shape
        # Assert
        status = False
        try:
            test = 'cfm_linear test'
            np.testing.assert_array_almost_equal(cfm_linear, testData['cfm_linear'], decimal=5)
            test = 'H test'
            np.testing.assert_array_almost_equal(H, testData['H'], decimal=5)
            test = 'b test'
            np.testing.assert_array_almost_equal(b, testData['b'], decimal=5)
            status = True
        except AssertionError as e:
            print e.message + '\n\nFailed at %s in test_get_linear_constraints' %(test)
            pass
        self.assertTrue(status)

    def test_bundletoa_update(self):
        # Fixture
        testData = sio.loadmat(os.path.join(testDataDir,'data_bundletoa_update'))
        D = testData['D']
        I = testData['I'] - 1 # Different indexation in Python and Matlab
        J = testData['J'] - 1 # Different indexation in Python and Matlab
        xt = testData['xt']
        yt = testData['yt']
        
        # Fixture
        res, jac = util.calcresandjac(D, I, J, xt, yt)
        dz_A = -(jac.T.dot(jac) + np.eye(jac.shape[1]))
        dz_b = jac.T.dot(res)
        
        case = 'multiply_LS'
        if case == 'multiply_LS':
            dz = scipy.linalg.solve(dz_A, dz_b)
        if case == 'sparse_LS':
            rcond = 1e-15
            dz = scipy.sparse.linalg.lsqr(dz_A, dz_b, rcond)[0]
            dz = np.reshape(dz, (len(dz),1))
        if case == 'numpy_LS':
            dz = numpy.linalg.lstsq(dz_A, dz_b)[0]
        if case == 'sparse_mult':
            A = scipy.sparse.csr_matrix(dz_A)
            b = scipy.sparse.csr_matrix(dz_b)
            dz = scipy.sparse.linalg.spsolve(A, b)
            dz = np.reshape(dz, (len(dz),1))
        
        xtn, ytn = util.updatexy(xt, yt, dz)
        res2, jac2 = util.calcresandjac(D, I, J, xtn, ytn)

        # Assert            
        status = False
        test = None
        try:
            test = 'A formation test'
            np.testing.assert_array_almost_equal(dz_A, testData['dz_A'], decimal=5)
            test = 'B formation test'
            np.testing.assert_array_almost_equal(dz_b, testData['dz_b'], decimal=5)
            test = 'LS solution test'
            np.testing.assert_array_almost_equal(dz, testData['dz'], decimal=5)
            test = 'Update test'
            np.testing.assert_array_almost_equal(xtn, testData['xtn'], decimal=5)   
            np.testing.assert_array_almost_equal(ytn, testData['ytn'], decimal=5)   
            test = 'Updates jacobian test'
            np.testing.assert_array_almost_equal(jac2.toarray(), testData['jac2'], decimal=5) 
            test = 'Updates residual test'
            np.testing.assert_array_almost_equal(res2, testData['res2'], decimal=5) 
            status = True
        except AssertionError as e:
            print e.message + '\n\nFailed at %s in test_bundletoa_update' %(test)
            pass
        self.assertTrue(status) 

    def test_bundletoa(self):
        # Fixture
        testData = sio.loadmat(os.path.join(testDataDir,'data_bundletoa'))
        D = testData['D']
        I = testData['I'] - 1 # Different indexation in Python and Matlab
        J = testData['J'] - 1 # Different indexation in Python and Matlab
        x = testData['x']
        y = testData['y']
        print D.shape
        print J.shape
        print I.shape

        settings = Settings()
        settings.bundle.numberOfIterations = 30
        settings.bundle.counterLimit = 50
        settings.bundle.numericalLimit = 1e-4

        # Test
        xopt, yopt, res, jac = util.bundletoa(D, I, J, x, y, settings)

        # Assert            
        status = False
        test = None
        try:
            test = 'residual test'
            np.testing.assert_array_almost_equal(res, testData['res'], decimal=0)
            test = 'x optimal test'
            np.testing.assert_array_almost_equal(xopt, testData['xopt'], decimal=0)
            test = 'y optimal test'
            np.testing.assert_array_almost_equal(yopt, testData['yopt'], decimal=0)            
            status = True
        except AssertionError as e:
            print e.message + '\n\nFailed at %s test_bundletoa' %(test)
            pass
        self.assertTrue(status)

    def test_toa_3d_bundle(self):
        testData = sio.loadmat(os.path.join(testDataDir,'data_toa_3d_bundle'))
        d = testData['d']
        rin = testData['rin']
        sin = testData['sin']
        inliers = testData['inliers']

        settings = Settings()
        settings.bundle.numberOfIterations = 30
        settings.bundle.counterLimit = 50
        settings.bundle.numericalLimit = 1e-4

        # Assert
        rout, sout, _, _ = util.toa_3d_bundle(d, rin, sin, inliers, settings)

        # Test
        # Assert            
        status = False
        test = None
        try:
            test = 'r optimal test'
            np.testing.assert_array_almost_equal(rout, testData['rout'], decimal=1)
            test = 's optimal test'
            np.testing.assert_array_almost_equal(sout, testData['sout'], decimal=0)            
            status = True
        except AssertionError as e:
            print e.message + '\n\nFailed at %s test_bundletoa' %(test)
            pass
        self.assertTrue(status)

if __name__ == '__main__':
    # Fixture
    testData = sio.loadmat(os.path.join(testDataDir,'data_toa_normalize'))
    rin = testData['rin']
    sin = testData['sin']

    # Assert
    rout, sout = util.toa_normalise(rin, sin)
        
    unittest.main()