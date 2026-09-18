#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 3 08:43:46 2026
DepositData4

A class to study ellipsometry data on LiF deposits.
See the Ellipsometry1 JupyterLab notebook.

Version 2 can handle the more complex data in the new files.
The lines now contain

x
y
Film thickness (t1)
Film grading???? (grading)
Film roughness (roughness)
Interlayer thickness (t2)
Interlayer void % ??? (voidness)

NOTE that the x and y arrays start at the top right and scan left and
down! Actually it is worse than that, the x values zig-zag back and
forth across the space. I do make the assumption that the first point
will be an outer edge point.

Version 3 stores all data in complete maps rectangular arrays using masks to
mark data points that arre missing either systematically (edges) or randomly.
This allows operations between two data sets using the
intersection of the masks to guarantee only using fully valid points.

This means that building a DepositData is considerably more complex as
__init__ now has to scan through the input arrays copying valid data
and building the mask. I have also made the t2, grading, roughness, and
voidness arguments optional since they will not make sense for composite
data sets.

Version4 supports the creation of a data copy rotated +90°.
This means that we need two different ways of building a new DepositData.
If we are creating a new data set from information read from a file then
we are passed 3 or more 1-D arrays and we have to sort them into our arrays.
If we are making one DepositData as a (possibly modified) copy of another
then we are passed a set of pre-sorted arrays and we just copy them into place.
Unfortunately, both kinds of data are stored in 1-D numpy arrays so that
we can't tell the apart by inspection. Thus I have added a 'sort' flag to
tell the constructor how to handle the data. I have moved the sorting
algorithm into a separate function to make __init__ more readable.

@author: bcollett
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

class DepositData4:
    # Normal input is arrays of x, y, t, etc. points that lie on parts of a
    # square grid and have to be sorted into place. If sort is set to False
    # then the input arrays are just copied into place.
    #
    def __init__(self, x, y, t1=None, t2=None, 
                 grading=None, roughness=None, voidness=None, 
                 name=None, mask=None, sort=True):
        self.name = name
        self.debugging = False
        #
        # Start by making sure that all the input data are the same
        # shape (if present).
        #
        good_shape = x.shape
        if y.shape != good_shape:
            raise ValueError('DepositData: y data shape must match x data shape.')
        if t1 is not None and t1.shape != good_shape:
            raise ValueError('DepositData: t1 data shape must match x data shape.')
        if t2 is not None and t2.shape != good_shape:
            raise ValueError('DepositData: t2 data shape must match x data shape.')
        if grading is not None and grading.shape != good_shape:
            raise ValueError('DepositData: grading data shape must match x data shape.')
        if roughness is not None and roughness.shape != good_shape:
            raise ValueError('DepositData: roughness data shape must match x data shape.')
        if voidness is not None and voidness.shape != good_shape:
            raise ValueError('DepositData: voidness data shape must match x data shape.')
        #
        # Create the data structure to be filled later.
        # This makes sure that all the members get created.
        #
        self.mask = None
        self.xdata = None
        self.ydata = None
        self.t1 = None
        self.t2 = None
        self.grading = None
        self.roughness = None
        self.voidness = None
        self.zdata = None
        self.fit_z = None
        self.resids = None
        #
        # Do we sort the data?
        #
        if sort:
            # We have raw data that we must sort into 2-D arrays before copying
            self._sort_from(x, y, t1, t2, 
                            grading, roughness, voidness)
        else:
            # we have sorted data that we just need to copy
            if mask is None:
                raise ValueError('DepositData: mask is required if sort is False.')
            self._copy_from(mask, x, y, t1, t2, 
                            grading, roughness, voidness)
        
        #
        # Useful statistics
        #
        if self.zdata is not None:
            self.avg_thickness = np.average(self.zdata)
        else:
            self.avg_thickness = 0.0
    
    def set_name(self, new_name):
        self.name = new_name
    
    # new_rot90 makes a copy of the data rotated 90 degrees in the positive
    # (counter-clockwise) direction. It should rotate whatever data the original
    # posesses but initially rotates only the zdata. This means that each
    # x,y is moved to -y,x. Update to rotate all data.
    def new_rot90(self):
        result = DepositData4(self.xdata,
                              self.ydata,
                              self.t1,
                              self.t2,
                              self.grading,
                              self.roughness,
                              self.voidness,
                              mask=self.mask,
                              sort=False,
                              name=f'Rotated {self.name}')
        npoint = len(self.xdata)
        xmax = np.max(result.xdata)
        ymax = np.max(result.ydata)
        nside = int(np.sqrt(npoint))
        eps = 1.0e-4
        for i in range(npoint):
            newx = -result.ydata[i]
            newy = result.xdata[i]
            ny = nside - 1 - int((ymax - newy) / self.dy + eps)
            nx = nside - 1 - int((xmax - newx) / self.dx + eps) # Assumes square
#            print(f'{y[i]} -> {ny}; {x[i]} -> {nx}')
            index = ny * nside + nx
            result.mask[index] = self.mask[i]
            result.t1[index] = self.t1[i]
            result.zdata[index] = self.zdata[i]
            if self.t2 is not None:
                result.t2[index] = self.t2[i]
            if self.grading is not None:
                result.grading[index] = self.grading[i]
            if self.roughness is not None:
                result.roughness[index] = self.roughness[i]
            if self.voidness is not None:
                result.voidness[index] = self.voidness[i]
        return result
    
    # new_with_mask generates a new DepositData keeping only those data
    # points that satisfy a mask. The mask must be an np.array of logic
    # values the same size as the x,y, z, and n arrays. All entries where
    # the mask is True will kept and the rest discarded.
    # Since the data are already masked, this operation changes only
    # the mask of the data. All the rest are simply copied.
    def new_with_mask(self, mask):
        self.mask = self. mask & mask
        result = DepositData4(self.xdata,
                              self.ydata,
                              self.t1,
                              self.t2,
                              self.grading,
                              self.roughness,
                              self.voidness,
                              mask=self.mask,
                              sort=False)
        if self.name is not None:
            result.set_name(f'Masked{self.name}')
        return result
        
    # Plotting functions
    def Ron_plot_on(self, ax):
        c = ax.tricontourf(self.xdata[self.mask], self.ydata[self.mask], 
                           self.zdata[self.mask], levels=201, cmap='hsv')
        ax.scatter(self.xdata[self.mask], self.ydata[self.mask], 
                   c='k', marker='.', s=10)
        plt.colorbar(c, label='Thickness (nm)')
        plt.grid(True, linestyle='-', linewidth=1, alpha=0.7)
        plt.xlim(-4, 4)
        plt.ylim(-4, 4)
        plt.xlabel('X (cm)')
        plt.ylabel('Y (cm)')
        plt.title(f'LiF Thickness for {self.name}')

    def Ron_plot(self):
        ax = self.new_plot()
        self.Ron_plot_on(ax)
    
    def new_plot(self):
        fig, ax = plt.subplots(figsize=(8,6))
        return ax
    
    def resid_plot_on(self, ax):
        if self.resids is not None:
            c = ax.tricontourf(self.xdata, self.ydata, self.resids, levels=201, cmap='hsv')
#            ax.scatter(self.xdata, self.ydata, c='k', marker='.', s=10)
            plt.colorbar(c, label='Thickness (nm)')
            plt.grid(True, linestyle='-', linewidth=1, alpha=0.7)
            plt.xlim(-4, 4)
            plt.ylim(-4, 4)
            plt.xlabel('X (cm)')
            plt.ylabel('Y (cm)')
            plt.title(f'Residuals for {self.name}')

    def resid_plot(self):
        ax = self.new_plot()
        self.resid_plot_on(ax)

    def resid_plot3D_on(self, ax):
        if self.resids is not None:
            c = ax.tricontourf(self.xdata[self.mask], self.ydata[self.mask], self.resids[self.mask], c= self.resids[self.mask], levels=201, cmap='hsv')
#            ax.scatter(self.xdata, self.ydata, c='k', marker='.', s=10)
            plt.colorbar(c, label='Thickness (nm)')
            plt.grid(True, linestyle='-', linewidth=1, alpha=0.7)
            plt.xlim(-4, 4)
            plt.ylim(-4, 4)
            plt.xlabel('X (cm)')
            plt.ylabel('Y (cm)')
            plt.title(f'Residuals for {self.name}')

    def resid_plot3D(self):
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surf = ax.scatter(self.xdata, self.ydata, self.resids, c= self.resids, cmap='hsv')

    # Contour plot an array with same size as our data
    def contour_plot(self, array):
        fig, ax = plt.subplots(figsize=(8,6))
#        c = ax.tricontourf(self.xdata, self.ydata, array, levels=201, cmap='hsv')
        c = ax.tricontourf(self.xdata[self.mask], self.ydata[self.mask], 
                           array[self.mask], levels=201, cmap='hsv')
        plt.colorbar(c, label='Thickness (nm)')
        plt.grid(True, linestyle='-', linewidth=1, alpha=0.7)
        plt.xlim(-4, 4)
        plt.ylim(-4, 4)
        plt.xlabel('X (cm)')
        plt.ylabel('Y (cm)')
        

    # Fit to some functional form
    # The fitting function should have the signature
    # fit_func(params, arg) where <params> will be the current
    # values of the fitting parameters and <arg> a tuple containing
    # the class. The function should return a single float value to minimize.
    # The p0 argument should hold the initial values of the parameters.
    def fit_to(self, fit_func, p0):
        result = optimize.minimize(fit_func, p0, args=(self, ), method='BFGS', options={'maxiter': 5000})
        fit_func(result.x, self)
        return result
    
    # Functions to evaluate various shapes given a parameter array.
    # These can be bound with fit_helper to give specific fitting functions.
    #
    # First one is an angled plane.
    # parameter definitions
    # p[0]: z0
    # p[1]: theta in radians
    # p[2]: the slope
    def angled_plane(self, p):
        stheta = np.sin(p[1])
        ctheta = np.cos(p[1])
        zplane = p[0] - p[2] * (self.xdata*ctheta + self.ydata*stheta)
        return zplane

    # Next an angled quadratic.
    # parameter definitions
    # p[0]: z0
    # p[1]: theta in radians
    # p[2]: the slope
    # p[3]: the quadraticity
    def angled_hump(self, p):
        stheta = np.sin(p[1])
        ctheta = np.cos(p[1])
        xi = self.xdata*ctheta + self.ydata*stheta
        zplane = p[0] - p[2] * xi - p[3] * xi * xi
        return zplane


    # Next an angled pure quartic.
    # parameter definitions
    # p[0]: z0
    # p[1]: theta in radians
    # p[2]: the slope
    # p[3]: the quarticity
    def angled_ridge(self, p):
        stheta = np.sin(p[1])
        ctheta = np.cos(p[1])
        xi = self.xdata*ctheta + self.ydata*stheta
        zplane = p[0] - p[2] * xi - p[3] * (xi ** 4)
        return zplane

    # Then an angled even quartic.
    # parameter definitions
    # p[0]: z0
    # p[1]: theta in radians
    # p[2]: the slope
    # p[3]: the quarticity
    # p[4]: the quadraticity
    def angled_quartic(self, p):
        stheta = np.sin(p[1])
        ctheta = np.cos(p[1])
        xi = self.xdata*ctheta + self.ydata*stheta
        zplane = p[0] - p[2] * xi - p[3] * (xi ** 2)- p[4] * (xi ** 4)
        return zplane
    
    #
    # Private functions used as helpers by e.g. the constructor
    #
    # _sort_from builds new arrays from unsorted data. It is basically
    # the old constructor
    #
    def _sort_from(self, x, y, t1, t2=None, 
                 grading=None, roughness=None, voidness=None):
        #
        # Extract grid parameters from x and y arrays.
        # NOTE: At least at the moment I am going to assume a square
        # grid with dx = dy
        # Use these params to build x and y arrays.
        #
#        print(f'x[0] - x[1] = {x[0]} - {x[1]}  = {x[0] - x[1]}')
        self.dx = x[0] - x[1]
        self.dy = self.dx
        self.radius = np.sqrt(y[0] * y[0] + x[0] * x[0] + 0.1)
        ymax = y[0]
#        print(f'ymax = {ymax}')
        
        self.xdata, self.ydata = DepositData4.new_grid(self.radius, self.dx)
#        print(self.xdata.shape)
        #
        # Alloc. all storage.
        # First required.
        self.t1 = np.zeros_like(self.xdata)
        self.zdata = np.zeros_like(self.xdata)
        self.mask = np.zeros_like(self.xdata, dtype=bool)
        nside = int(np.sqrt(self.xdata.shape[0]))
#        print(f'nside = {nside}')
        #
        # Then the optional ones
        # Note constructor already set all to None 
        if t1 is not None:
#            print('Using t1')
            self.t1 = np.zeros_like(self.xdata)
        else:
            self.t2 = None
        if t2 is not None:
#            print('Using t2')
            self.t2 = np.zeros_like(self.xdata)
        else:
            self.t2 = None
        #
        if grading is not None:
            self.grading = np.zeros_like(self.xdata)
        else:
            self.grading = grading
        #
        if roughness is not None:
            self.roughness = np.zeros_like(self.xdata)
        else:
            self.roughness = roughness
        #
        if voidness is not None:
            self.voidness = np.zeros_like(self.xdata)
        else:
            self.voidness = voidness

        #
        # Copy valid elements into arrays and build mask.
        #
        mask_count = 0;
        eps = 1.0e-4
        good_shape = x.shape
        for i in range(good_shape[0]):
#        for i in range(15):
#            print(f'{y[i]}, {ymax - y[i]}, {int((ymax - y[i]) / self.dy)}')
            ny = nside - 1 - int((ymax - y[i]) / self.dy + eps)
            nx = nside - 1 - int((ymax - x[i]) / self.dx + eps) # Assumes square
#            print(f'{y[i]} -> {ny}; {x[i]} -> {nx}')
            index = ny * nside + nx
            self.mask[index] = True
            if np.abs(x[i] - self.xdata[index]) > eps:
                print(f'x error at i = {i}, x={x[i]}, index={index}, xdata={self.xdata[index]}')
                continue
            if np.abs(y[i] - self.ydata[index]) > eps:
                print(f'y error at i = {i}, y={y[i]}, index={index}, ydata={self.ydata[index]}')
                continue
            mask_count += 1
            #
            # Point is valid. Copy data from i to index
            #
            if t1 is not None:
                self.t1[index] = t1[i]
                self.zdata[index] = self.t1[i]
                if t2 is not None:
                    self.t2[index] = t2[i]
                    self.zdata[index] = t2[i] + t1[i]
            if grading is not None:
                self.grading[index] = grading[i]
            if roughness is not None:
                self.roughness[index] = roughness[i]
            if grading is not None:
                self.roughness[index] = roughness[i]
            if index < 0: 
                print(f'{i} -> {index}')
                print(f'i = {i}, x,y={x[i]}, {y[i]}, '
                      f'index={index}, xdata,ydata={self.xdata[index]},{self.ydata[index]};'
                      f' z={self.zdata[index]}')
#                print(f'i = {i}, y={y[i]}, index={index}, ydata={self.ydata[index]}')
    
    #
    # _copy_from builds new arrays as copies of existing ones.
    #
    def _copy_from(self, mask, x, y, t1, t2=None, 
                 grading=None, roughness=None, voidness=None, zdata=None):
        self.mask = mask.copy()
        self.xdata = x.copy()
        self.ydata = y.copy()
        if t1 is not None:
            self.t1 = t1.copy()
        if t2 is not None:
            self.t2 = t2.copy()
        if grading is not None:
            self.grading = grading.copy()
        if roughness is not None:
            self.roughness = roughness.copy()
        if voidness is not None:
            self.voidness = voidness.copy()
        if zdata is not None:
            self.zdata = zdata.copy()
        else:
            if t1 is not None:
                if t2 is not None:
                    self.zdata = t1 + t2
                else:
                    self.zdata = t1.copy()
            else:
                self.zdata = np.zeros_like(self.xdata)

    # Class functions
    # Class factory function to build a DepositData from a file
    def from_file(filename, data_name = None):
        arr = np.loadtxt(filename, delimiter='\t', skiprows=5)
        dd = DepositData4(arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 5], arr[:, 3], arr[:, 4], arr[:, 6]) 
        if data_name is not None:
            dd.set_name(data_name)
        return dd
    
    # helper for fits. Bind to a real method name before use.
    def fit_helper(params, args, m_name):
        data = args
        func = getattr(data, m_name)
        data.fit_z = func(params)
        data.resids = data.zdata - data.fit_z
        res_sq = data.resids[data.mask] ** 2
        return np.sum(res_sq)

    # Helper to build the complete x and y arrays required to make the
    # Data complete. It includes all points that are integer multiples
    # of dx or dy, symmetrically about zero, within +/-radius in each direction.
    def new_grid(radius, dx, dy = None):
        if dy is None:
            dy = dx
#        rad_sq = float(radius) ** 2
#        print(f'new_grid: dx = {dx}, rad = {radius}')
#        print('rad_sq =', rad_sq)
            
        num_xsegs = int(radius / dx)
        nXSeg = 2 * num_xsegs + 1
        newx = np.zeros(nXSeg * nXSeg)

        num_ysegs = int(radius / dy)
        nYSeg = 2 * num_ysegs + 1
        newy = np.zeros(nYSeg * nYSeg)
#        print(f'nx = {nXSeg}, ny = {nYSeg}')
#        print(newx.shape, newy.shape)
        index = 0
 #       deleted = 0
        for j in range(nYSeg):
            y = (j - num_ysegs) * dy
            for i in range(nXSeg):
                x = (i - num_xsegs) * dx
#                rsq = x * x + y * y
#                if rsq <= rad_sq:
                newx[index] = x
                newy[index] = y
                index += 1
#                else:
#                    deleted += 1
#                    if deleted < 10:
#                        print(x, y)
#        print(f'Deleted {deleted}')
#        print(newx[:30])
#        print(newy[:30])
        return np.round(newx[:index],4), np.round(newy[:index],4)


if __name__ == '__main__':
    s3 = DepositData4.from_file('Run6_26_Scan3.txt', 'Scan3')
    s3.Ron_plot()
    r3 = s3.new_rot90()
    r3.Ron_plot()
    