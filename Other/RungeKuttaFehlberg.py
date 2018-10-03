#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 20:26:05 2018

@author: rivertz
"""

import numpy as np
import math as m 
import sys
import timeit
import matplotlib.pyplot as plt

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
                 stepsize, tolerance):
        self.F = function;
        self.dim = dimension;
        self.h = stepsize;
        self.tol = tolerance;
    
    def step(self, W_in):
        s = np.zeros((6, self.dim));

        for i in range(0, 6):
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

def F(Y):
    M = np.array([[1, 1],
                [-1, 1]]);
    res = np.ones(3);
    res[1:3] = M.dot(Y[1:3]);
    return res;
def F_test(Y):
    M = np.array([[0.49119653, 0.32513304, 0.98057799],
                [0.20768544, 0.97699416, 0.18220559],
                [0.96407071, 0.18373237, 0.95307793]]);
    res = np.ones(4);
    res[1:4] = M.dot(Y[1:4]);
    return res;

def y_1(t):
    return m.e**t * m.cos(t);

def y_2(t):
    return -m.e**t * m.sin(t);

def main():
    """ 6.3.1.a
      h = 0.25, [0,1]
      y'1 = y_1 + y_2
      y'2 = −y_1 + y_2
      y_1(0) = 1
      y_2(0) = 0

      y_1(t) = e**t * cos(t), y_2(t) = -e**t * sin(t)
      """

    W = np.array([0, 1, 0]);
    E = 0;
    h = 0.25;
    tol = 05e-14;
    tEnd = 1.0;

    rkf54 = RungeKuttaFehlberg54(F, 3, h, tol);
    while (W[0] < tEnd):
        W, E = rkf54.safeStep(W);

    rkf54.setStepLength(tEnd - W[0]);
    W, E = rkf54.step(W);

    print(W, E);
    print("Total error: ", [abs(W[1] - y_1(1)), abs(W[2] - y_2(1))]);

def tid(tol):
    W = np.array([0, 1, 1, 1]);
    h = 0.1;
    tEnd = 2.0;
    rkf54 = RungeKuttaFehlberg54(F_test, 4, h, tol);

    while (W[0] < tEnd):
        W, E = rkf54.safeStep(W);

    rkf54.setStepLength(tEnd - W[0]);
    W, E = rkf54.step(W);
    return W, E

def tidplot():
    ekspo = np.arange(0, 20, 1)
    tiden = np.ones(20)

    for i in range(20):
        tol = 1 * 10 ** -i
        tiden[i] = timeit.timeit('tid({})'.format(tol), 'from __main__ import tid', number=10) / 10;
        # print('Toleranse: {}\t Tid: {}'.format(tol,tiden[i]));

    plt.plot(ekspo, tiden, 'b-')
    plt.ylabel('Tid per rkf45 (s)')
    plt.xlabel('Toleransens eksponent (10^-x)')
    plt.show()
    
if __name__ == "__main__":
    # execute only if run as a script
    main();
    tidplot();



