#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 20:26:05 2018

@author: rivertz
"""

import math as m
import sys
import time

import numpy as np
import matplotlib.pyplot as plt


class RungeKuttaFehlberg54:
    A = np.array(
        [[0,     0,     0,    0,   0, 0],
         [1/4,     0,     0,    0,   0, 0],
         [3/32,     9/32,     0,    0,   0, 0],
         [1932/2197, -7200/2197,  7296/2197,    0,   0, 0],
         [439/216,    -8,  3680/513, -845/4104,   0, 0],
         [-8/27,     2, -3544/2565, 1859/4104, -11/40, 0]])
  
    B = np.array(
        [[25/216,    0,   1408/2565,  2197/4104, -1/5, 0],
         [16/135,    0, 6656/12825, 28561/56430, -9/50, 2/55]])

    def __init__(self, f, dimension,
                 step_size, tolerance):
        self.f = f
        self.dim = dimension
        self.h = step_size
        self.tol = tolerance
    
    def step(self, w_in):
        s = np.zeros((6, self.dim))

        for i in range(0, 6):
            s[i, :] = self.f(w_in + self.h * self.A[i, 0:i].dot(s[0:i, :]))

        z_out = w_in + self.h * (self.B[0, :].dot(s))
        w_out = w_in + self.h * (self.B[1, :].dot(s))

        error = np.linalg.norm(w_out - z_out, 2) / np.linalg.norm(w_out, 2)
        return w_out, error

    def safe_step(self, w_in):
        w_out, error = self.step(w_in)
        # Check if the error is tolerable
        if not self.is_error_tolerated(error):
            # Try to adjust the optimal step length
            self.adjust_step(error)
            w_out, error = self.step(w_in)
        # If the error is still not tolerable
        counter = 0
        while not self.is_error_tolerated(error):
            # Try if dividing the steplength with 2 helps. 
            self.divide_step_by_two()
            w_out, error = self.step(w_in)
            counter += 1
            if counter > 10:
                print("System is unreliable, terminating.")
                sys.exit(-1)
            
        self.adjust_step(error)
        
        return w_out, error

    def is_error_tolerated(self, error):
        return error < self.tol

    def adjust_step(self, error):
        if error == 0:
            s = 2
        else:
            s = m.pow(self.tol * self.h / (2 * error), 0.25)
        self.h = s * self.h

    def divide_step_by_two(self):
        self.h = self.h / 2
        
    def set_step_length(self, step_length):
        self.h = step_length        
 
        
"""def F(Y):
    M = np.array([[0.49119653, 0.32513304, 0.98057799],
                [0.20768544, 0.97699416, 0.18220559],
                [0.96407071, 0.18373237, 0.95307793]])
    res = np.ones(4)
    res[0:3] = M.dot(Y[0:3])
    return res"""


def func(y):
    matrix = np.array([[1, 1], [-1, 1]])
    res = np.ones(3)
    res[1:3] = matrix.dot(y[1:3])
    return res


def y_1(t):
    return m.e**t * m.cos(t)


def y_2(t):
    return -m.e**t * m.sin(t)


def main(tolerance):
    """ 6.3.1.a
    h = 0.25, [0,1]
    y'1 = y_1 + y_2
    y'2 = −y_1 + y_2
    y_1(0) = 1
    y_2(0) = 0

    y_1(t) = e**t * cos(t), y_2(t) = -e**t * sin(t)
    """
    w = np.array([0, 1, 0])
    h = 0.25
    tol = tolerance
    t_end = 1.0
    rkf54 = RungeKuttaFehlberg54(func, 3, h, tol)

    while w[0] < t_end:
        w, error = rkf54.safe_step(w)
        
    rkf54.set_step_length(t_end - w[0])
    w, error = rkf54.step(w)

    result = [w[1], w[2]]
    total_error = [abs(w[1] - y_1(1)), abs(w[2] - y_2(1))]
    return result, total_error
    
    
if __name__ == "__main__":
    # execute only if run as a script
    data = [
        [],  # Tolerance
        [],  # Runtime
        [],  # Result
        [],  # Total error
    ]

    for i in range(0, 20):
        tolerance = 10**-i
        t1 = time.time()
        result, total_error = main(tolerance)
        t2 = time.time()
        runtime = t2-t1
        data[0].append(tolerance)
        data[1].append(runtime)
        data[2].append(result)
        data[3].append(total_error)

    fig, ax1 = plt.subplots()
    ax1.plot(data[0], data[1], "b.", lw=0, ms=6)
    ax1.set_xlabel("Tolerance")
    ax1.set_ylabel("Runtime (s)", color="b")
    ax1.tick_params("y", colors="b")

    ax2 = ax1.twinx()
    ax2.plot(data[0], data[3], "r.", lw=0, ms=3)
    ax2.set_ylabel("Error", color="r")
    ax2.tick_params("y", colors="r")

    fig.tight_layout()

    plt.xscale("log")
    plt.grid(True)
    plt.show()
