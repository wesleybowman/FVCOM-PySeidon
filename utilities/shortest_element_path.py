from __future__ import division
import numpy as np
import scipy.spatial
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
import seaborn

class shortest_element_path:
    '''
Description:
----------
A class/structure using FVCOM data to construct a path across the centers of
the grid provided by the FVCOM class.

Inputs:
------

ind = [-66.3419, -66.3324, 44.2755, 44.2815]

test = FVCOM(filename, ax=ind)

path = shortest_element_path(test.latc,test.lonc,test.lat,test.lon,test.nv,test.h)

elements, coordinates = path.getTargets([[41420,39763],[48484,53441],
                                        [27241,24226],[21706,17458]])

path.graphGrid(plot=True)


Options:
-------
debug = True will print out a lot of useful information.

Notes:
-----

'''
    def __init__(self, lonc, latc, lon, lat, nv, h, debug=False):

        self.debug = debug
        #self.data = nc.Dataset(filename,'r')

        #latc = self.data.variables['latc'][:]
        #lonc = self.data.variables['lonc'][:]
        self.lonc = lonc
        self.latc = latc
        self.lat = lat
        self.lon = lon
        self.nv = nv
        self.h = h

        z = np.vstack((lonc,latc)).T

        self.points = map(tuple,z)

        print 'File Loaded'

        # make a Delaunay triangulation of the point data
        self.delTri = scipy.spatial.Delaunay(self.points)
        print 'Delaunay Triangulation Done'

        # create a set for edges that are indexes of the points
        self.edges = []
        self.weight = []
        # for each Delaunay triangle
        for n in xrange(self.delTri.nsimplex):
            # for each edge of the triangle
            # sort the vertices
            # (sorting avoids duplicated edges being added to the set)
            # and add to the edges set

            self.edge = sorted([self.delTri.vertices[n,0], self.delTri.vertices[n,1]])
            a = self.points[self.edge[0]]
            b = self.points[self.edge[1]]
            self.weight = (np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2))
            self.edges.append((self.edge[0], self.edge[1],{'weight':self.weight}))

            self.edge = sorted([self.delTri.vertices[n,0], self.delTri.vertices[n,2]])
            a = self.points[self.edge[0]]
            b = self.points[self.edge[1]]
            self.weight = (np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2))
            self.edges.append((self.edge[0], self.edge[1],{'weight':self.weight}))


            self.edge = sorted([self.delTri.vertices[n,1], self.delTri.vertices[n,2]])
            a = self.points[self.edge[0]]
            b = self.points[self.edge[1]]
            self.weight = (np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2))
            self.edges.append((self.edge[0], self.edge[1],{'weight':self.weight}))

        print 'Edges and Weighting Done'

        # make a graph based on the Delaunay triangulation edges
        self.graph = nx.Graph(self.edges)
        #print(graph.edges())

        print 'Graph Constructed'

        self.pointIDXY = dict(zip(range(len(self.points)), self.points))

    def getTargets(self, source_target, coords=False):
        ''' getTargets takes in a source_target. This can be in the form of
        actual coordinates or indexes that correspond to coordinates. If
        coords=True, then the source_target needs to be in EXACT coordinates.
        The suggestion method would be to run closest_point on the coordinates
        that are wanted, and then pass the indexes into source_target'''

        self.elements = []
        self.coordinates = []
        self.maxcoordinates = []
        self.mincoordinates = []
        for i in source_target:
            source = i[0]
            target = i[1]

            if self.debug:
                print '\n'
                print 'Source'
                print source

            s = source
#
            if self.debug:
                print 'Target'
                print target
            t = target

            if coords:
                for key, value in self.pointIDXY.items():
                    if value==source:
                        print 'Source'
                        print key
                        s = key

                    if value==target:
                        print 'Target'
                        print key
                        t = key

            #print s,t
            shortest = nx.shortest_path(self.graph,source=s,target=t,weight='weight')
#            dist = nx.shortest_path_length(self.graph,source=s,target=t,weight='weight')

            if self.debug:
                print 'Shortest Path (by elements)'
                print shortest

            self.elements.append(shortest)

            coords = [self.pointIDXY[i] for i in shortest]
            self.coordinates.append(coords)
            self.maxcoordinates.append(np.max(np.array(coords), axis=0))
            self.mincoordinates.append(np.min(np.array(coords), axis=0))

#            print 'Shortest Distance (by coordinates)'
#            print dist

        return self.elements, self.coordinates

    def graphGrid(self,narrowGrid=False, plot=False):
        ''' A method to graph the grid with the shortest path plotted on the
        grid. narrowGrid will limit the field of view down to only show the
        paths drawn and not the entire grid. If only one path is drawn, this
        can skew the way the grid looks, and is sometime better to view the
        whole grid and zoom in. The plot option is in cause you want to choose
        when to plot the graph in a different part of the code.'''
        lat = self.lat
        lon = self.lon
        nv = self.nv.T - 1
        h = self.h

        tri = Tri.Triangulation(lon, lat, triangles=nv)  # xy or latlon based on how you are #Grand Passage

        levels=np.arange(-38,6,1)   # depth contours to plot

        fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')
        ax = fig.add_subplot(111,aspect=(1.0/np.cos(np.mean(lat)*np.pi/180.0)))
        plt.tricontourf(tri, -h,levels=levels,shading='faceted',cmap=plt.cm.gist_earth)
        plt.triplot(tri)
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.gca().patch.set_facecolor('0.5')
        cbar=plt.colorbar()
        cbar.set_label('Water Depth (m)', rotation=-90,labelpad=30)

        scale = 1
        ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        ax.xaxis.set_major_formatter(ticks)
        ax.yaxis.set_major_formatter(ticks)
        plt.grid()

        maxlat, maxlon = np.max(self.maxcoordinates,axis=0)
        minlat, minlon = np.min(self.mincoordinates,axis=0)
        if narrowGrid:
            ax.set_xlim(minlon,maxlon)
            ax.set_ylim(minlat,maxlat)


        zz = len(self.elements)
        for i,v in enumerate(self.elements):
            source = self.pointIDXY[v[0]]
            target = self.pointIDXY[v[-1]]
            lab = '({:.6},{:.6})-({:.6},{:.6})'.format(source[0], source[1],
                                                       target[0], target[1])

            plt.scatter(self.lonc[v], self.latc[v],
                        s=80, label=lab, c=plt.cm.Set1(i/zz))

        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=2, ncol=3)
        if plot:
            plt.show()


if __name__ == '__main__':

    filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'

    test = shortest_element_path(filename)

    test.getTargets([[41420,39763],[48484,53441],[27241,24226],[21706,17458],[14587,5416]])
    test.graphGrid(narrowGrid=True)

    element_path, coordinates_path = test.getTargets([[41420,39763]])
