from __future__ import print_function, absolute_import
import numpy as np
from scipy.spatial import Delaunay 
import scipy.optimize as op
import functools 
import logging
import string

import gravpy.trinterior as trint
import gravpy.plots

class gravlens:
        
    def __init__(self,carargs,polargs,modelargs,show_plot=True,include_caustics=True,image=np.random.normal(0,0.25,2),recurse_depth=3,caustics_depth=8,logging_level='info'):
        self.carargs = carargs
        self.xspacing = carargs[0][2]
        self.yspacing = carargs[1][2]
        self.polargs = polargs
        self.modelargs = modelargs
        self.show_plot = show_plot
        self.include_caustics = include_caustics
        self.image = image
        self.recurse_depth = recurse_depth
        self.caustics_depth= caustics_depth
        self.cache = {}

        self.validate_arguments()

        
        levels = {'critical':50, 'error':40, 'warning':30, 'info':20, 'debug':10, 'notset':0}
        if isinstance(logging_level,str):
            logging_level = levels[logging_level.lower()]
        else:
            assert isinstance(logging_level,int)
            
        logging.basicConfig(level=logging_level,format='%(message)s')
        self.log = logging.getLogger(__name__)
        
        self.num_eval = 0
        self.log.info("")

    def validate_arguments(self):
        def is_a_sequence_type(arg):
            return isinstance(arg, (list, tuple, np.ndarray))

        
        if not np.array(self.carargs).shape==(2,3):
            raise AssertionError("'carargs' are not of the correct shape")
        if not is_a_sequence_type(self.polargs[0]):
            raise AssertionError("'polargs' is not a list of list(s)")
        if not is_a_sequence_type(self.modelargs):
            raise AssertionError("'modelargs' is not a list of models")
        
    
        
    def relation(self,x,y):
        '''tells us if the point pair (x,y) is outside, inside, or on the critical curve'''
        
        dif = self.magnification(x,y)
        return np.sign(dif) # -1 for inside, +1 for outside, 0 for exactly on

    @staticmethod
    def polartocar(r,th):
        '''convert polar coordinates to cartesian coordinates'''
        r  = np.array(r)
        th = np.array(th)
        
        return np.transpose(r*[np.cos(th),np.sin(th)])
    
    def mapping(self,v):
        '''Function used for root finding. \'v\' is the true position in the image plane (what we want to find). Not intended for vectorized inputs. Returns the deflection vector along with the inverse-magnification matrix (aka the Jacobian).'''
    
        w,z = self.image
        x,y = v #actual position we want to find
    
        phiarr = self.potdefmag(x,y,numexpr=False).ravel()
        phix,phiy,phixx,phiyy,phixy = phiarr[1:6]
    
        return [[x-w-phix,y-z-phiy],# deflection
                [[1-phixx, -phixy], # jacobian
                 [-phixy ,1-phiyy]]]
    
    
    def carmapping(self,x,y):
        '''mapping of cartesian coordinates from image to source plane'''
        self.log.debug("******Mapping Call******")
        phiarr = self.potdefmag(x,y)
        phix,phiy = phiarr[1:3]
        
        return np.transpose([x-phix,y-phiy])
    
    def magnification(self,x,y):
        '''returns the magnification of the points '''
        self.log.debug("***Magnification Call***")
        phiarr = self.potdefmag(x,y)
        phixx,phiyy,phixy = phiarr[3:6]
        
        return (1-phixx)*(1-phiyy)-phixy**2

    def memoize(func):
        
        @functools.wraps(func)
        def memoizer(*args, **kwargs):
            __self,xi,yi = args
            cache = __self.cache
            
            #turn scalars into vectors of length 1
            x = np.atleast_1d(xi) 
            y = np.atleast_1d(yi) 

            xc = []
            yc = []
            use_cache = set() 
            for i in xrange(x.size):
                key = make_key(x[i],y[i])
                if key in cache:
                    use_cache.add(i)
                else:
                    xc.append(x[i])
                    yc.append(y[i])
                        
                
            __self.log.info("{:>13} in cache, calculating {:<8} new value(s)"
                               .format(str(len(x)-len(xc)) + "/" + str(len(x)),len(xc)))

            if len(xc) is not 0:
                phiarray = np.transpose(func(__self,xc,yc,**kwargs))

            
            outarray = []
            counter = iter(range(len(xc)))
            for i in xrange(x.size):
                key = make_key(x[i],y[i])
                if i in use_cache:
                    outarray.append(cache[key])
                else:
                    arr = phiarray[counter.next()]
                    outarray.append(arr)
                    cache[key] = arr
            
            return np.transpose(outarray)

        def make_key(x,y):
            """Define the key string for the cache dictionary"""
            return str(x) + " , " + str(y)
        
        return memoizer
    
    @memoize
    def potdefmag(self,xi,yi,numexpr=True):
        '''The wrapper used to find the phi values given models' parameters. The output is (6,x) where x is the length of the x,y arguments given in the invocation. This command seeks out the correct module to contact for each model calculation.'''
        phi2Darray = []
    
        #turn scalars into vectors of length 1
        x = np.atleast_1d(xi) 
        y = np.atleast_1d(yi)

        # make sure both are same size
        
        if not x.shape == y.shape:
            raise ValueError("x and y arrays are not the same size.")

        self.num_eval += x.size

        self.log.debug("Evaluating %d point(s)..." % x.size)

        for mass_component in self.modelargs:
            
            phiarray = mass_component.phiarray(x,y,numexpr=numexpr) 
            phi2Darray.append(phiarray)
                
        return np.nan_to_num(np.sum(phi2Darray,axis=0))
    
    
    def mag_of_cells(self,cells):
        '''Takes a list of cells and returns the magnification values for each point in the cells. Retains shape and order of the original list of cells.'''
        num_points_in_cell = cells.shape[1] #usually 4, unless we're caching mag values (then it's 3)
        cells_x_y = cells.reshape((-1,2))
        cells_x = cells_x_y[:,0]
        cells_y = cells_x_y[:,1]
    
        mag_x_y = self.relation(cells_x,cells_y)
        
        return np.reshape(mag_x_y,(-1,num_points_in_cell))
    
    @staticmethod
    def cell_mag_change(cells_mag):
        '''Takes a list of magnification values of a list of cells and returns a boolean mask for which the magnification changes across a cell. Uses numpy vectorization.'''
        fir,sec,thr,frt = cells_mag.T

        with np.errstate(invalid='ignore'): 
            less1 = np.vstack((fir*sec,sec*frt,frt*thr,thr*fir)) < 1
        
        output = np.any(less1,axis=0)
                
        if np.count_nonzero(output) == 0: #np.all doesn't catch the error
            raise ValueError("Magnification does not change across grid, change grid parameters so that critical curves are seen.")
        
        else:
            return output

    
    def subdivide_cells(self,cells,cell_depth):
        '''Divides each cell in a list of cells (a cell given by [p1, p2, p3, p4]) into four smaller cells. On the cartesian grid, the points are ordered as so:
        p2--------p4         |
        |         |      II  |  I
        |         |     -----|-----
        |         |      III |  IV
        p1--------p3         |
        '''
        #cry/reimplement if we need to subdivide cells by any different than 4:1

        dx = self.xspacing/ 2**cell_depth
        dy = self.yspacing/ 2**cell_depth
        
        #*** The following code is not used, but helps illustrate what the output is ***
        # below code uses broadcasting to 'shrink' each cell into a fourth of its size, but still retaining one of its vertices. This happens four times, shrinking each time towards one of the four vertices, leaving four quarter cells that together make up the original cell.
        # quadrant1 = cells + [[dx, dy],[dx,  0],[0  , dy],[0  ,  0]]
        # quadrant2 = cells + [[0 , dy],[0 ,  0],[-dx, dy],[-dx,  0]]
        # quadrant3 = cells + [[0 ,  0],[0 ,-dy],[-dx,  0],[-dx,-dy]]
        # quadrant4 = cells + [[dx,  0],[dx,-dy],[0  ,  0],[0  ,-dy]]

        """
        The following code computes the intermediate points (p0,p24,p12,p34,p13) and then combines them into the four quadrants.
        p2---p24--p4         |
        |    |    |      II  |  I
        p12--p0--p34    -----|-----
        |    |    |      III |  IV
        p1--p13---p3         |
        """
        
        p1,p2,p3,p4 = np.transpose(cells,[1,0,2])
        p24 = p4 + [-dx,  0]
        p12 = p2 + [0  ,-dy]
        p13 = p1 + [dx ,  0]
        p34 = p3 + [0  , dy]
        p0  = p1 + [ dx, dy]

        quadrant1 = [p0,p24,p34,p4]
        quadrant2 = [p12,p2,p0,p24]
        quadrant3 = [p1,p12,p13,p0]
        quadrant4 = [p13,p0,p3,p34]

        # 2nd element: need x,y values for computing mag values
        return (np.transpose(np.hstack((quadrant1,quadrant2,quadrant3,quadrant4)), [1,0,2]),
                np.transpose(np.vstack((p0,p12,p13,p34,p24)))) 
    
    def for_points5_wrapper_cached(self,cells,mag_cells,cell_depth):
        '''Function that subdivdes given cells, computes the magnification values for new points, merges the magnification values from \'mag_cells\' with the newly computed magnifcation values, and returns the magnification values and cells where a critical curve was detected.'''
                
        subdivided_cells, points_lists = self.subdivide_cells(cells,cell_depth)
        
        m0,m12,m13,m34,m24 = np.split(self.relation(points_lists[0],points_lists[1]),5)
        m1,m2,m3,m4 = mag_cells.T # old mag values (did not change)

        # look in subdivide_cells() for this recipe
        m_quadrant1 = [m0,m24,m34,m4]
        m_quadrant2 = [m12,m2,m0,m24]
        m_quadrant3 = [m1,m12,m13,m0]
        m_quadrant4 = [m13,m0,m3,m34]

        mag_combined = np.transpose(np.hstack((m_quadrant1,m_quadrant2,m_quadrant3,m_quadrant4)))

        # which cells had a change in magnification
        mag_change_mask = self.cell_mag_change(mag_combined) 
        
        # return mag values of selected cells (used for the next level in gridding) and the cells themselves that had a mag change
        return [np.compress(mag_change_mask,mag_combined    ,axis=0),
                np.compress(mag_change_mask,subdivided_cells,axis=0)]
    
    def points5(self,xran,yran):
        '''A vectorized approach to bulding a 'recursive' subgrid without recursion. Algorithm works by vectorizing each level of cell-size, handling each level in one complete calculation before proceeding to the next. '''
                
        x = xran[0:-1]
        y = yran[0:-1]
        xs, ys = np.tile(x,len(y)),np.repeat(y,len(x))
    
        grid_pairs = np.column_stack((np.tile(xran,len(yran)),np.repeat(yran,len(xran))))
    
        gridm_x_y = np.vstack((np.dstack((xs,xs,xs+self.xspacing,xs+self.xspacing)),np.dstack((ys,ys+self.yspacing,ys,ys+self.yspacing))))
        cells = np.transpose(gridm_x_y,[1,2,0])
    
        # we don't want to subdivide the first iteration
        cells_mag = self.mag_of_cells(cells)
        mag_change_mask = self.cell_mag_change(cells_mag)

        #equivalent to cells[mag_change_mask] but faster
        cells_sel = np.compress(mag_change_mask,cells,axis=0) 

        # ""; = cells_mag[mag_change_mask]
        cells_mag = np.compress(mag_change_mask,cells_mag,axis=0) 
        output_pairs = []
        
        depth = self.caustics_depth if self.include_caustics else self.recurse_depth
        for i in range(depth):
            
            cells_mag,cells_sel = self.for_points5_wrapper_cached(cells_sel,cells_mag,i+1)
                        
            output_pairs.append(cells_sel)
    
        
        output_pairs = np.vstack(output_pairs).reshape((-1,2))
        if self.include_caustics:
            # don't want the vertices of each cell; just the (center) of each cell
            self.critical_lines = np.transpose(np.mean(cells_sel,axis=1)) 
            
        return np.vstack((grid_pairs,output_pairs))
        
    
    
    def generate_ranges(self):
        '''Generates the sequences used for the other core & gridding functions.
        Returns an array containing:
        [ [x,y], [r,theta], [critx, crity], [causticsx, causticsy]]
        where critx,crity refers to the x and y values for the critical curves,
        causticsx, causticsy refers to the x and y values for the caustics,
        x,y are the cartesian ranges for the originial mesh-grid,
        r,theta are the cartesian x and y values for the supplementary polar grid(s).
        If caustics is set to False, then [critx, crity] and [causticsx, causticsy] are not returned-- Instead just [x,y]  and [r,theta] are returned.
        '''
        
        # initial cartesian grid, coarse,
        [xlowerend, xupperend, xspacing],[ylowerend,yupperend,yspacing] = self.carargs
        
        self.x = np.arange(xlowerend,xupperend+xspacing,xspacing)
        self.y = np.arange(ylowerend,yupperend+yspacing,yspacing)
    
        #supplemental polar grid(s)
        polargrid = [] 
        
        for grid in self.polargs:
            center, rupper, rdivisions, thetadivisions = grid

            delzero = +1
            logr1 = np.log10(rupper+delzero)
    
            r = np.logspace(0,logr1,num=rdivisions)-delzero
            theta = np.linspace(0,2*np.pi,thetadivisions,endpoint=False)
            
            ## the rs and thetas formed from a cartesian product. Used for vector operations.
            temprs = np.tile(r,len(theta))
            tempthetas = np.repeat(theta,len(r))
    
            shiftedpairs = self.polartocar(temprs,tempthetas) + center
    
            polargrid.append(shiftedpairs)
        
        self.polargrid = np.vstack(polargrid)
    
    
    def transformations(self):
        '''Generates the subgridding points (more points around the critical curve), the transformations from image plane to source plane, and the Delaunay Triangulization object for plotting.'''
        
        # generate subgrid on cartesian grid
        cargrid = self.points5(self.x,self.y)
        
        # combine list of cartesian and polar pairs
        self.grid = np.concatenate((cargrid,self.polargrid),axis=0)

        # transform pairs from image to source plane
        self.transformed = self.carmapping(self.grid[:,0],self.grid[:,1]) 
        
        if self.include_caustics:
            criticalx, criticaly = self.critical_lines
            causticsx, causticsy = np.transpose(self.carmapping(criticalx,criticaly))
            self.caustics = [[criticalx,criticaly],[causticsx,causticsy]]
        else:
            self.caustics = None
                
        
    def find_source(self):
        '''Employs the algorithm in the 'trinterior' module to find the positions of the image in the image plane. Returns the coordinate pair(s) in an array.'''
        
        # generate Delaunay object for triangulization/triangles
        dpoints = Delaunay(self.grid) 
        self.simplices = dpoints.simplices
    
        lenstri = np.take(self.transformed,self.simplices,axis=0) #triangles on the source plane
        imagetri= np.take(self.grid       ,self.simplices,axis=0) #triangles on the image plane

        try:
            # list of which triangles contain point on source plane
            indices = trint.find2(self.image,lenstri) 
        except ValueError:
            self.log.fatal("Image location corresponding to source position not found. Re-run with bigger grid or/and use a different source position.")
            self.realpos = []
            return
        
        sourcetri = imagetri[indices]
        
        # list of the centroid coordinates for the triangles which contain the point 'image'
        sourcepos = np.mean(sourcetri,axis=1)
        
        # use centroid coordinates as guesses for the actual root finding algorithm
        realpos = [(op.root(self.mapping,v,jac=True,tol=1e-4)).x
                   for v in sourcepos]
        
        
        self.realpos = np.array(realpos)
    
    def run(self):
        '''The master command that wraps and executes all the commands to run a gridding example. Use this function (excusively) when using this module.'''
        
        self.generate_ranges()
            
        self.transformations()
        
        self.find_source()
        
        if self.show_plot:
            self.plot()

        self.log.info("%d points evaluated" % self.num_eval)
        self.log.info("%d points in cache " % len(self.cache.keys()))

        
    def plot(self):
 
        [xlowerend, xupperend, xspacing],[ylowerend,yupperend,yspacing] = self.carargs
        
        plots.source_image_planes(
                self.grid,self.transformed,self.simplices,
                self.realpos,self.image,
                xlowerend,xupperend,ylowerend,yupperend,
                caustics=self.caustics)
    
    
    
    
