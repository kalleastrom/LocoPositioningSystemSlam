"""
A module for automatic anchor localization in an UWB system
"""
import numpy as np
from numpy.linalg import cholesky, inv, eig, solve, det, svd
from numpy.matlib import repmat
import os, json
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D # Only used for 3D plotting

class AnchorLocalizer(object):
    """
    The anchor localizer object takes a data file or direct measurements in an
    UWB positioning system and computes the position of the anchors based
    on ranging measurements between the UAV and the UWB anchors.
    
    The object can be set to include a UAV model to improve the anchor position
    estimate with the known system dynamics. However, the data used to
    propagate the dynamics in time must then be time stamped and logged at the
    same instance in time as the ranging measurement.
    
    """
    def __init__(self, dataDirectory):
        """
        Constructor of the anchor localizer object
        """
        self.dataDirectory           = dataDirectory
        self.numberOfAnchors         = None
        self.numberOfMeasurements    = None
        self.information             = None
        self.ranges                  = None
        self.initialPositions        = None
        self.finalPositions          = None
        self.finalPositionDeviations = None

        # Intermediary attributes for debugging 
        self.D = 3  # The number of dimensions

    def __str__(self):
        """
        String representation of the anchor localizer object
        
        :returns: representation
        :rtype: string
        """
        return "Anchor localizer object"

    def __call__(self):
        """
        Calls the solver and finds the anchor positions
        """
        # Parameters
        self.numIterInit = 2 # Number of iterations to find an initial solution
        
        # Finds an initial feasible solution from the ranging data
        self.find_feasible_solution()

        # Optimization to estimate the anchor positions
        r, s, res, jacobian = self.find_optimal_solution()

    def load(self, filename):
        """
        Loads a data file as a JSON object

        The data must be on the form of a dictionary with the keys

            {
                'ranges':...,
                'acceleration':...,
                'rates':...,
                'thrust':...,
                'torque':...,
                'times':...
             }

        The 'ranges' data must be a sequence of M measured distances from N 
        anchors, stored in a list of lists with size (N,M), such that each row
        contins all ranging information from a single anchor. This field is
        mandatory.

        All other fields are optional and can be included as required. If used,
        these fields must be filled with measurements taken at the same time
        as the ranging is logged, with the time stamp fields detailing the time
        at which measurements were taken. 

        See /examples for further information.

        :param filename: Name of the file to load
        :type filename: string
        """
        path = os.path.join(self.dataDirectory, filename)
        try:
            with open(path) as trajectoryFile:
                data = json.load(trajectoryFile)
        except:
            raise Exception(('Could not load the data file "%s" from the'+
                             'directory %s' % (filename)))
        self.ranges = np.array(data['ranges'])
        self.numberOfAnchors = self.ranges.shape[0]
        self.numberOfMeasurements = self.ranges.shape[1]
        self.information = np.array(data['info'])

    def verbose(self):
        """
        Prints the statistics of the most recent anchor localization
        
        :param state: The state of the object
        :type state: string
        """
        state = ''
        if self.information is None:
            state += '* No trajectory loaded'
        else:
            state += '* Loaded: %s\n' % (self.information) 
            if self.ranges is None:
                state += '* No syntacticallty correct ranging data found\n'
            else:
                state += '* Ranging data loaded of dimensions %s\n' % str(self.ranges.shape)
            if self.initialPositions is None:
                state += '* No solution has been computed\n'
            else:
                state += '* A solution has been found\n' % str(self.ranges.shape)
        return state

    def plot(self, options):
        """
        Visualizes the anchor localization in space

        The method is called with an options argument which is a list
        containing the requested plots, where
        * "ranges" - Visualises the measured ranges between the anchors and the
            UAV in a 2D image with distance as a function of the anchor and
            the measurements
        * "anchors" - Visualises anchors with initial guess, an estimated
            position and the associated standard deviation of each estimate.

        :param options: A list with string valued options
        :type options: [string]
        """
        if 'ranges' in options:
            if not self.ranges is None:
                fig = plt.figure(1)
                plt.matshow(self.ranges,
                            cmap='gray',
                            extent=[0,2*self.numberOfAnchors,0.5,self.numberOfAnchors+0.5])
                plt.title('The measured range [m] of an anchor at a measurement i')
                plt.xlabel('Masurements')
                plt.ylabel('Anchor')
                plt.colorbar()
                plt.show()
        plt.cla()
        if 'anchors' in options:
            if self.initialPositions is None or self.finalPositionDeviations is None:
                self.message('No initial anchor position has been computed')
            else:
                fig = plt.figure(5)
                ax = fig.gca(projection='3d')
                ax.plot(self.initialPositions[0,:],
                        self.initialPositions[1,:],
                        self.initialPositions[2,:], 'gx')
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_zlabel('z')
                for ii in range(self.finalPositions.shape[1]):
                    ax.plot(self.finalPositions[0,:],
                            self.finalPositions[1,:],
                            self.finalPositions[2,:], 'r+')
                    x, y, z = self._get_ellipsoid(center = self.finalPositions[:,ii],
                                                  radius = self.finalPositionDeviations[:,ii])
                    ax.plot_surface(x, y, z,
                                    rstride=4,
                                    cstride=4,
                                    color='b',
                                    linewidth=0,
                                    antialiased=False,
                                    alpha=0.1)
                plt.show()
    
    def message(self, string):
        """
        A method that is called when communicating with the user
        
        :returns: message
        :rtype: string
        """
        print string
    
    def _get_ellipsoid(self, center, radius):
        """
        Computes an N-point ellipsoid with given radii and center for plotting
        
        :returns: x, y, z
        :rtype: np.array(N), np.array(N), np.array(N)
        """
        N = 50
        u = np.linspace(0, 2 * np.pi, N)
        v = np.linspace(0, np.pi, N)
        x = center[0] + radius[0] * np.outer(np.cos(u), np.sin(v))
        y = center[1] + radius[1] * np.outer(np.sin(u), np.sin(v))
        z = center[2] + radius[2] * np.outer(np.ones_like(u), np.cos(v))
        return x, y, z
    
    ###########################################################################
    # Functions used to localize the anchors, refer to the publication below
    # and use the bib-tex entry for citation
    #
    # @inbook{
    #     6505724206a7482c8bec0010f5a11197,
    #     title = "Robust Time-of-Arrival Self Calibration with Missing Data and Outliers",
    #     author = "Batstone, {Kenneth John} and Magnus Oskarsson and Karl \AA str\o" m",
    #     year = "2016",
    #     month = "9",
    #     doi = "10.1109/EUSIPCO.2016.7760673",
    #     pages = "2370--2374",
    #     booktitle = "2016 24th European Signal Processing Conference (EUSIPCO)",
    #     publisher = "Institute of Electrical and Electronics Engineers Inc.",
    #     address = "United States"
    # }
    ###########################################################################
    
    def find_feasible_solution(self):
        """
        Find an initial feasible solution to the problem using RANSAC bundling
        """
        self.ransac_5_rows();
        self.bundle_rank();
        for kk in range(self.numIterInit):
            self.ransac_more_rows()
            self.bundle_rank()
            self.ransac_more_cols()
            self.bundle_rank()

    def find_optimal_solution(self):
        """
        Find an initial feasible solution to the problem using L2 optimization
        """
        # Go from the found Bhat to {H,b} to {r,s}
        u, s, v = svd(self.Bhat[1:,1:])
        s = np.diag(s)
        xr  = np.transpose(u[:,0:self.D])
        yr  = np.dot(s[0:self.D,0:self.D], np.transpose(v[:,0:self.D]))
        xtp = np.concatenate((np.zeros((self.D,1)), xr),1)
        yt  = np.concatenate((np.zeros((self.D,1)), yr),1)
        xt  = xtp/-2.0;
        
        # Get the linear constraints
        H, b = self.get_linear_constraints(xt, yt)
        
        # Guarantee positive definiteness of H in order to cholesky factorize
        realEig = np.real(eig(H)[0]) # The real part of the eigenvalues of H
        if min(realEig) > 0:
            L = cholesky(inv(H))
        else:
            H = H + (-min(realEig) + 0.1) * np.eye(3)
            L = cholesky(inv(H))

        # Finds the R and S vectors
        r00 = solve(np.transpose(L), xt)
        s00 = np.dot(L, (yt + repmat(b, 1, self.Bhat.shape[1])));
        
        # Reshapes the initial r00 and s00 vectors TODO: Why is this done?
        # TODO: write out equations
        r0 = r00
        s0 = s00
        
        # Finds the bundled solution and normalizes
        r1, s1, res1, jacobian1 = self.bundle(r0, s0)
        r2, s2 = self.normalize(r1, s1)
        
        # Finds the smoothened solution
        r, s, res, jacobian = self.bundle_with_smoother(r2, s2)

        # computes the covariance matrix and checks that it exists
        #J = jacobian[:,0:2]
        #JT = np.transpose(jacobian[:,0:2])
        #JTJ = np.dot(JT, J)
        #if abs(det(JTJ)) < 1.0e-6:
        #    self.issue_warning("The covariance matrix is singular to working precision")
        # TODO Compute covariance matrix C
            
        if self.debug:
            # Outputs all temporary variables when debugging or testing
            return (H, b, xr, yr, xtp, yt, xt,
                    r00, s00, r1, s1, res1, jacobian1,
                    r2, s2, r, s, res, jacobian)
        return r, s, res, jacobian
        
    def get_linear_constraints(self, xt, yt):
        """
        Sets up the linear constraints (not using the multipol package)
        """
        Bhatcol1 = self.Bhat[:,0]
        #todo write out equations
        H = np.zeros((3,3))
        b = np.zeros((3,1))
        return H, b

        
    def ransac_5_rows(self):
        """
        
        """
        print 'Not yet implemented'

    def ransac_more_rows(self):
        """
        
        """
        print 'Not yet implemented'

    def ransac_more_cols(self):
        """
        
        """
        print 'Not yet implemented'

    def bundle_rank(self):
        """
        
        """
        print 'Not yet implemented'


    def bundle(self, r, s):
        """
        
        """
        print 'Not yet implemented'
        r, s, res, jacobian = 0,0,0,0
        return r, s, res, jacobian

    def normalize(self, r, s):
        """
        
        """
        print 'Not yet implemented'
        r, s  = 0,0
        return r, s

    def bundle_with_smoother(self, r, s):
        """
        
        """
        print 'Not yet implemented'
        r, s, res, jacobian = 0,0,0,0
        return r, s, res, jacobian

class DiscreteDynamics(object):
    
    def __init__(self, Nx, Nu, Ny, Nc):
        """
        Contructor of the discrete dynamics super class
        
        :param Nx: The number of states (x)
        :param Nu: The number of control signals (u)
        :param Ny: The number of measurement signals (y)
        :param Nc: The number of constant terms (c)
        :type Nx: int
        :type Nu: int
        :type Ny: int
        :type Nc: int
        """
        self.Nx = Nx                           # Number of states
        self.Nu = Nu                           # Number of control signals
        self.Ny = Ny                           # Number of measurement signals
        self.Nc = Nc                           # Number of constant signals
        self.Nexpm = 3                         # Number of terms in the disc.
        self.Ad = np.zeros((self.Nx,self.Nx))  # Discrete time state matrix
        self.Bd = np.zeros((self.Nx,self.Nu))  # Discrete time c.sig. matrix
        self.Cd = np.zeros((self.Ny,self.Nx))  # Discrete time measurement matrix
        self.x = np.zeros((self.Nx,1))         # State vector
        self.u = np.zeros((self.Nu,1))         # Control signal vector
        self.y = np.zeros((self.Ny,1))         # Control signal vector

    def __call__(self, u, h):
        """
        Propagates the discrete time model forward by a time step h
        
        :param u: Control signal vector at the current time t = hk
        :param h: The time step in seconds
        :type u: np.array((Nc,1))
        :type h: int
        """
        self.u = u
        self.x = np.dot(self.Ad, self.x) + np.dot(self.Bd, self.u) + self.Gd
        self.y = np.dot(self.Cd, self.x)
        return self.x, self.y

    def discretize_system(self, h):
        """
        Discretises a continuous time model at a time step of h with ZOH
        
        :param h: The time step in seconds
        :type h: int
        """
        dim = self.Nx + self.Nu + self.Nc
        
        # Forms a concatenated matrix on the form tM = [A, B, G; 0, 0, 0]
        totalMatrix = np.zeros((dim,dim))
        totalMatrix[0:self.Nx, 0:self.Nx] = self._get_A()
        totalMatrix[0:self.Nx, self.Nx:(self.Nx+self.Nu)] = self._get_B()
        totalMatrix[0:self.Nx, (self.Nx+self.Nu):dim] = self._get_G()
        
        # Discretises system using the exponential matrix formulation
        expM = np.eye(dim)
        expMdt = np.eye(dim)
        for ii in range(1, self.Nexpm):
            expMdt = np.dot((h * totalMatrix / ii), expMdt)
            expM += expMdt

        # Extract discrete time system matrices
        self.Ad = expM[0:self.Nx, 0:self.Nx]
        self.Bd = expM[0:self.Nx, self.Nx:(self.Nx+self.Nu)]
        self.Gd = expM[0:self.Nx, (self.Nx+self.Nu):dim]
    
    def set_state(self, x):
        """ Sets the system state vector externally """
        self.x = x

    def set_control_signals(self, u):
        """ Sets the system control signal vector externally """
        self.u = u

    def _get_A(self):
        """
        Defined in inheriting classes to return the continuous time
        system matrix, could be made a non-linear of the current x and u
        """
        raise Exception('_get_A should be overshadwed in inheriting objects')
        A = np.zeros((self.Nx,self.Nu))
        return A
    
    def _get_B(self):
        """
        Defined in inheriting classes to return the continuous time control
        signal matrix, could be made a non-linear of the current x and u
        """
        raise Exception('_get_B should be overshadwed in inheriting objects')
        B = np.zeros((self.Nx,self.Nu))
        return B

    def _get_G(self):
        """
        Defined in inheriting classes to return the continuous time
        constant matrix, could be made a non-linear of the current x and u
        """
        raise Exception('_get_G should be overshadwed in inheriting objects')
        G = np.zeros((self.Nx,self.Nc))
        return G

    def _get_C(self):
        """
        Defined in inheriting classes to return the continuous time
        measurement matrix, could be made a non-linear of the current x and u
        """
        raise Exception('_get_C should be overshadwed in inheriting objects')
        C = 0.0
        return C

    def __str__(self):
        return 'A discrete dynamics object'

class TrippleIntegrator(DiscreteDynamics):
    
    def __init__(self, dt):
        super(TrippleIntegrator, self).__init__(Nx=9, Nu=3, Ny=3, Nc=1)
        self.discretize_system(dt)

    def _get_A(self):
        A = np.zeros((self.Nx,self.Nx))
        A[0:6, 3:9] = np.eye(6)
        return A
    
    def _get_B(self):
        B = np.zeros((self.Nx,self.Nu))
        B[6:9,0:3] = np.eye(3)
        print B
        return B

    def _get_G(self):
        G = np.zeros((self.Nx,self.Nc))
        return G

    def __str__(self):
        return 'A three dimensional tripple integrator object'