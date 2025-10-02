#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 11:45:50 2024

@author: kendrickshepherd
"""

import sys
import numpy as np
import math
import matplotlib
from matplotlib import pyplot as plt

# IMPORT (or copy) your code from HW3 here
# which evaluated a Lagrane polynomial basis function
# import <UNCOMMENT THIS LINE AND INPUT YOUR FILENAME HERE>
sys.path.append('../HW3_Univariate_Lagrange/')

import Univariate_Lagrange_Basis_Functions as HW3

# higher-dimensional basis function with multi-index
def MultiDimensionalBasisFunctionIdxs(idxs,degs,interp_pts,xis):
    val = 1
    for i in range(0,len(idxs)):
        val *= HW3.LagrangeBasisEvaluation(degs[i],interp_pts[i],xis[i],idxs[i])
    return val
    
# given an integer and a vector of basis function degrees extract
# a corresponding multi-index out
def ExtractMultiIndexFromSingleIndex(A,degs):
    bfs = degs[0]+1
    horiz = A % (bfs)
    idxs = [horiz]
    for i in range(1,len(degs)):
        idxs.append(A // bfs)
        bfs *= (degs[i]+1)
    return idxs

# higher-dimensional basis function with single index
def MultiDimensionalBasisFunction(A,degs,interp_pts,xis):
    # bfs = degs[0]+1
    # horiz = A % (bfs)
    # idxs = [horiz]
    # for i in range(1,len(degs)):
    #     idxs.append(A // bfs)
    #     bfs *= (degs[i]+1)
    idxs = ExtractMultiIndexFromSingleIndex(A, degs)
    return MultiDimensionalBasisFunctionIdxs(idxs,degs,interp_pts,xis)
    
# plot of 2D basis functions with A a single index
def PlotTwoDimensionalParentBasisFunction(A,degs,npts = 101,contours = True):
    interp_pts = [np.linspace(-1,1,degs[i]+1) for i in range(0,len(degs))]
    xivals = np.linspace(-1,1,npts)
    etavals = np.linspace(-1,1,npts)
    
    Xi,Eta = np.meshgrid(xivals,etavals)
    Z = np.zeros(Xi.shape)
    
    for i in range(0,len(xivals)):
        for j in range(0,len(etavals)):
            Z[i,j] = MultiDimensionalBasisFunction(A, degs, interp_pts, [xivals[i],etavals[j]])
    
    # contour plot
    if contours:
        fig, ax = plt.subplots()
        surf = ax.contourf(Eta,Xi,Z,levels=100,cmap=matplotlib.cm.jet)
        ax.set_xlabel(r"$\xi$")
        ax.set_ylabel(r"$\eta$")
        fig.colorbar(surf)
        plt.show()
    # 3D surface plot
    else:
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surf = ax.plot_surface(Eta, Xi, Z, cmap=matplotlib.cm.jet,
                       linewidth=0, antialiased=False)
        ax.set_xlabel(r"$\xi$")
        ax.set_ylabel(r"$\eta$")
        ax.set_zlabel(r"$N(\xi,\eta)$")
        plt.show()
        