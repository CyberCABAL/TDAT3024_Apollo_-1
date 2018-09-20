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
        self.state = np.array([[t],
                               [[objects[i].position[j] for i in range(self.objLen)] for j in range(self.dim)],
                               [[objects[i].dirVec[j] for i in range(self.objLen)] for j in range(self.dim)],
                               ]);

    def addObject(self, obj):
        self.objects.append(obj);
        self.objLen += 1;

    def stateToObjects(self):
        for i in range(self.objLen):
            objects[i].position = self.state[1][i];
            objects[i].dirVec = self.state[2][i];

    def velocityVector(self, i = None, name = None):
        if (i is None):
            if (name is not None):
                for j in range(self.objLen):
                    if (self.objects[j].name == name):
                        i = j;
            else:
                return np.zeros(3);
            
        dist = [(fy.dist(self.objects[i], self.objects[j]) if (i != j) else 0.1) for j in range(self.objLen)];      
        p0 = self.objects[i].position;

        vec = np.zeros(3);
        for j in range(self.objLen):
            if (i != j):
                r = dist[j];
                if (r == 0):
                    r = 0.1;
                pj = self.objects[j];
                vec += (pj.gm * (pj.position - p0)) / (r * r * r);
        return vec;

    def updateState(self, state):
        self.state = state;
        
    def ydot(self, x):
        return np.array([np.array([1]), self.state[1], [self.velocityVector(i) for i in range(self.objLen)]]);

    def step(self, rate):
        for i in range(rate - 1):
            self.estimate.safeStep(self.state);
        return self.estimate.safeStep(self.state);
                          
