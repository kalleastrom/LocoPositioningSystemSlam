"""
A module for automatic anchor localization in an UWB system
"""
import numpy as np
import json
import os
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
        # TODO: Find an initial feasible solution using the ransac algorithm

        # TODO: Iterate until a good solution is found

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
        self.numberOfMeasurements = self.ranges.shape[0]
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
                plt.matshow(self.ranges, cmap='gray',extent=[0,2,0,1])
                plt.title('The measured range [m] of an anchor at a measurement i')
                plt.xlabel('Masurements')
                plt.ylabel('Anchor')
                plt.colorbar()
                plt.show()
        if 'anchors' in options:
            if self.initialPositions is None or self.finalPositionDeviations is None:
                self.message('No initial anchor position has been computed')
            else:
                fig = plt.figure(2)
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