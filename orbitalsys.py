import fy
from RungeKuttaFehlberg import RungeKuttaFehlberg54

import numpy as np
import time
import matplotlib.pyplot as plot
import matplotlib.animation as animation

class System:
    def __init__(self,
                objects,
                 G = 6.67408 * 10**(-11), t = 0, estimate = None):
        self.GravConst = G;
        self.objects = objects;
        self.objLen = len(objects);
        self.dim = len(objects[0].position);
        self.state = np.array([np.array([t]),
                               np.array([np.array([objects[i].position[j] for i in range(self.objLen)]) for j in range(self.dim)]),
                               np.array([np.array([objects[i].dirVec[j] for i in range(self.objLen)]) for j in range(self.dim)]),
                               ]);
        self.estimate = estimate;

    def time_elapsed(self):
        return self.state[0][0];

    def objectsToState(self):
        self.state = np.array([self.state[0],
                               [[objects[i].position[j] for i in range(self.objLen)] for j in range(self.dim)],
                               [[objects[i].dirVec[j] for i in range(self.objLen)] for j in range(self.dim)],
                               ]);

    def addObject(self, obj):
        self.objects.append(obj);
        self.objLen += 1;

    def stateToObjects(self):
        state1 = self.state[1].T;
        state2 = self.state[2].T;
        for i in range(self.objLen):
            self.objects[i].position = state1[i];
            self.objects[i].dirVec = state2[i];

    def velocityVector(self, x, i = None, name = None):
        if (i is None):
            if (name is not None):
                for j in range(self.objLen):
                    if (self.objects[j].name == name):
                        i = j;
            else:
                return np.zeros(3);
        p = x[1];
        dist = self.distPos(p);
        p0 = p[:, i];

        vec = np.zeros(3);
        #print("p:", p);
        for j in range(self.objLen):
            if (i != j):
                r = dist[i][j];     # Can be simplified, there is only one i.
                if (r == 0):
                    r = 0.1;
                vec += (self.objects[j].gm * (p[:, j] - p0)) / (r * r * r);
        #print(vec);
        return vec;

    def distPos(self, x):
        return [[self.__tempSum(x, n, m)**(1.5) for m in range(self.objLen)] for n in range(self.objLen)];

    def __tempSum(self, x, n, m):
        sum0 = 0;
        for i in range(0, self.dim):
            temp = x[i][m] - x[i][n];
            sum0 += temp * temp;
        return sum0;

    def updateState(self, state):
        self.state = state;
        
    def ydot(self, x):
        v = np.array([self.velocityVector(x, i) for i in range(self.objLen)]);
        #print("Test1:", v);
        #print("Test2:", v.T);
        return np.array([np.array([1]), x[1], v.T]);

    #def ydot2(self, x, name):
    #    p = x[1];
    #    return np.array([np.array([1]), x[2], self.velocityVector(name)]);

    def step(self, rate = 3):
        W, E = self.estimate.safeStep(self.state);
        self.updateState(W);
        for i in range(rate - 1):
            W, E = self.estimate.safeStep(W);
            self.updateState(W);
        self.stateToObjects();
        return W, E;

    #def step2(self, rate):
    #    r = None;
    #    for j in range(objLen):
    #       w = np.array([self.state[0], self.state[1][j], self.state[2][j]]);
    #        for i in range(rate):
    #            r = self.estimate[j].safeStep(w);
    #    return r;
                          
