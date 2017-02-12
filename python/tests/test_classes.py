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
testDataDir = os.path.join(testDir, 'test-data')
sys.path.insert(0,modulesdDir)

import lps_utilities as util
from lps_classes import Settings, Solver, Solution

if sys.version_info < (3, 3):
    from mock import MagicMock
else:
    from unittest.mock import MagicMock

class TestSolver(unittest.TestCase):

    def setUp(self):
        self.sut = Solver()

    def test_that_solution_is_found(self):
        # Fixture
        
        #Assert
        
        #Test
        self.assertAlmostEqual(1,1)

class TestSolution(unittest.TestCase):

    def setUp(self):
        self.sut = Solution()

    def test_something(self):
        # Fixture
        
        #Assert
        
        #Test
        self.assertAlmostEqual(1,1)

class TestSettings(unittest.TestCase):

    def setUp(self):
        self.sut = Settings()

    def test_something(self):
        # Fixture
        
        #Assert
        
        #Test
        self.assertAlmostEqual(1,1)

if __name__ == '__main__':
    unittest.main()