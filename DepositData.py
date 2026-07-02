#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 09:54:46 2026

A class to study ellipsometry data on LiF deposits.
See the Ellipsometry1 JupyterLab notebook.

@author: bcollett
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

class DepositData:
    def __init__(self, x, y, z, n):
        self.name = None
        self.xdata = x
        self.ydata = y
        self.zdata = z
        self.ndata = n
        self.avg_thickness = np.average(self.zdata)
        self.fit_z = None
        self.resids = None
    
    def set_name(self, new_name):
        self.name = new_name
    
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
            c = ax.tricontourf(self.xdata, self.ydata, self.resids, levels=201, cmap='hsv')
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
        surf = ax.scatter(self.xdata, self.ydata, self.resids)
    
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

    # Class functions
    # Class factory function to build a DepositData from a file
    def from_file(filename, data_name = None):
        arr = np.loadtxt(filename, delimiter='\t', skiprows=5)
        dd = DepositData(arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]) 
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
