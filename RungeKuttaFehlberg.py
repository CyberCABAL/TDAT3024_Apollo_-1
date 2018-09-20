#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 20:26:05 2018

@author: rivertz
"""

import numpy as np
import math as m 
import sys

class RungeKuttaFehlberg54:
    A = np.array(
        [[   0     ,     0     ,     0      ,    0     ,   0   , 0],
         [   1/4   ,     0     ,     0      ,    0     ,   0   , 0],
         [   3/32  ,     9/32  ,     0      ,    0     ,   0   , 0],
         [1932/2197, -7200/2197,  7296/2197 ,    0     ,   0   , 0],
         [ 439/216 ,    -8     ,  3680/513  , -845/4104,   0   , 0],
         [  -8/27  ,     2     , -3544/2565 , 1859/4104, -11/40, 0]]);
  
    B = np.array(
        [[  25/216,    0     , 1408/2565 ,  2197/4104 , -1/5 , 0],
         [  16/135,    0     , 6656/12825, 28561/56430, -9/50, 2/55]]);

    def __init__(self, function, dimension,
                 stepsize, tolerance, n):
        self.F = function;
        self.dim = dimension;
        self.h = stepsize;
        self.tol = tolerance;
        self.n = n;
    
    def step(self, W_in):
        #s = np.zeros((6, self.dim));
        s = np.array([[np.array([0]), np.array([np.zeros(self.n) for i in range(self.dim)]), np.array([np.zeros(self.n) for i in range(self.dim)])] for j in range(6)]);

        for i in range(0, 6):
            print(W_in, "\n");
            print(self.A[i, 0:i], "\n");
            print(s[0:i, :], "\n");
            print(self.A[i, 0:i].dot(s[0:i, :]), "\n");
            print(W_in + self.h * self.A[i, 0:i].dot(s[0:i, :]));
            s[i, :] = self.F(W_in + self.h * self.A[i, 0:i].dot(s[0:i, :]));

        Z_out = W_in + self.h * (self.B[0, :].dot(s));
        W_out = W_in + self.h * (self.B[1, :].dot(s));

        E = np.linalg.norm(W_out - Z_out, 2) / np.linalg.norm(W_out, 2);
        return W_out, E;

    def safeStep(self, W_in):
        W_out, E = self.step(W_in);
        # Check if the error is tolerable
        if(not self.isErrorTolerated(E)):
            #Try to adjust the optimal step length
            self.adjustStep(E);
            W_out, E = self.step(W_in);
        # If the error is still not tolerable
        counter = 0;
        while(not self.isErrorTolerated(E)):
            #Try if dividing the steplength with 2 helps. 
            self.divideStepByTwo();
            W_out, E = self.step(W_in);
            counter += 1;
            if(counter > 10):
                print("System is unreliable, terminating.")
                sys.exit(-1);
            
        self.adjustStep(E);
        
        return W_out, E;

    def isErrorTolerated(self, E):
        return E < self.tol;

    def adjustStep(self, E):
        if(E == 0):
            s = 2;
        else:
            s = m.pow(self.tol * self.h / (2 * E), 0.25);
        self.h = s * self.h;

    def divideStepByTwo(self):
        self.h = self.h / 2;
        
    def setStepLength(self, stepLength):
        self.h = stepLength;

