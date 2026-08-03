#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 3 08:43:46 2026

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

@author: bcollett
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

class DepositData:
    def __init__(self, x, y, t1, t2, grading, roughness, voidness):
        self.name = None
        self.xdata = x
        self.ydata = y
        self.t1 = t1
        self.t2 = t2
        self.grading = grading
        self.roughness = roughness
        self.voidness = voidness
        self.zdata = self.t1+self.t2
        self.avg_thickness = np.average(self.zdata)
        self.fit_z = None
        self.resids = None
    
    def set_name(self, new_name):
        self.name = new_name
    
    # new_with_mask generates a new DepositData keeping only those data
    # points that satisfy a mask. The mask must be an np.array of logic
    # values the same size as the x,y, z, and n arrays. All entries where
    # the mask is True will kept and the rest discarded.
    def new_with_mask(self, mask):
        newx = self.xdata[mask]
        newy = self.ydata[mask]
        newz = self.zdata[mask]
        result = DepositData(newx, newy, newz, newn)
        if self.name is not None:
            result.set_name(f'Masked{self.name}')
        return result
        
    # Plotting functions
    def Ron_plot_on(self, ax):
        c = ax.tricontourf(self.xdata, self.ydata, self.zdata, levels=201, cmap='hsv')
        ax.scatter(self.xdata, self.ydata, c='k', marker='.', s=10)
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
            c = ax.tricontourf(self.xdata, self.ydata, self.resids, c= self.resids, levels=201, cmap='hsv')
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
        c = ax.tricontourf(self.xdata, self.ydata, array, levels=201, cmap='hsv')
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
    # p[3]: the quadraticity
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
    # p[3]: the quadraticity
    def angled_quartic(self, p):
        stheta = np.sin(p[1])
        ctheta = np.cos(p[1])
        xi = self.xdata*ctheta + self.ydata*stheta
        zplane = p[0] - p[2] * xi - p[3] * (xi ** 2)- p[4] * (xi ** 4)
        return zplane
    
    

    # Class functions
    # Class factory function to build a DepositData from a file
    def from_file(filename, data_name = None):
        arr = np.loadtxt(filename, delimiter='\t', skiprows=5)
        dd = DepositData(arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 5], arr[:, 3], arr[:, 4], arr[:, 6]) 
        if data_name is not None:
            dd.set_name(data_name)
        return dd
    
    # helper for fits. Bind to a real method name before use.
    def fit_helper(params, args, m_name):
        data = args
        func = getattr(data, m_name)
        data.fit_z = func(params)
        data.resids = data.zdata - data.fit_z
        res_sq = data.resids * data.resids
        return np.sum(res_sq)
