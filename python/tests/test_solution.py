import unittest
import scipy.io as sio
import os, sys, inspect
import numpy as np

curdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
rootdir = os.path.dirname(curdir)
modulesdDir = os.path.join(rootdir, 'modules')
dataDir = os.path.join(rootdir, 'data')
sys.path.insert(0,modulesdDir)

from locolib import AnchorLocalizer
testFile = 'TestData.mat'



class TestAnchorLocalizer(unittest.TestCase):

    def get_data_from_file(self, fn):
        try:
            data = sio.loadmat(fn)
        except:
            raise Exception(("Could not load file, check that %s exists.") % (fn))
        return data
    
    def test_find_optimal_solution(self):
        """ Tests the optimization against known data """
        al = AnchorLocalizer(dataDir)
        al.debug = True
        data = self.get_data_from_file(testFile)

        # Set up a known initial feaisble solution
        al.Bhat = data['Bhat']

        # Find anchors
        (H, b, xr, yr, xtp, yt, xt,
         r00, s00, r1, s1, res1, jacobian1,
         r2, s2, r, s, res, jacobian) = al.find_optimal_solution()

        # Verify solution agains Matlab data
        self.assertAlmostEqual(xr.all(),        data['xr'].all())
        self.assertAlmostEqual(yr.all(),        data['yr'].all())
        self.assertAlmostEqual(xtp.all(),       data['xtp'].all())
        self.assertAlmostEqual(xt.all(),        data['xt'].all())
        self.assertAlmostEqual(yt.all(),        data['yt'].all())
        self.assertAlmostEqual(H.all(),         data['H'].all())
        self.assertAlmostEqual(b.all(),         data['b'].all())
        self.assertAlmostEqual(r00.all(),       data['r00'].all())
        self.assertAlmostEqual(s00.all(),       data['s00'].all())
        self.assertAlmostEqual(r1.all(),        data['r1'].all())
        self.assertAlmostEqual(s1.all(),        data['s1'].all())
        self.assertAlmostEqual(res1.all(),      data['res1'].all())
        self.assertAlmostEqual(jacobian1.all(), data['jacobian1'].all())
        self.assertAlmostEqual(r2.all(),        data['r2'].all())
        self.assertAlmostEqual(s2.all(),        data['s2'].all())
        self.assertAlmostEqual(r.all(),         data['r'].all())
        self.assertAlmostEqual(s.all(),         data['s'].all())
        self.assertAlmostEqual(res.all(),       data['res'].all())
        self.assertAlmostEqual(jacobian.all(),  data['jacobian'].all())
 
    def test_strings_a_3(self):
        self.assertEqual('aaa', 'aaa')
 
if __name__ == '__main__':
    unittest.main()