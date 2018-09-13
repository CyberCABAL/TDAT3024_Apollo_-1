import fy
from RungeKuttaFehlberg import RungeKuttaFehlberg54

import numpy as np
import time
import matplotlib.pyplot as plot
import matplotlib.animation as animation

class System:
    def __init__(self,
                objects,
                 G = 6.67408 * 10**(-11), t = 0):
        self.GravConst = G;
        self.objects = objects;
        n = len(objects);
        dim = len(objects[0].position);
        self.state = np.array([[t],
                               [[objects[i].position[j] for i in range(n)] for j in range(dim)],
                               [[objects[i].dirVec[j] for i in range(n)] for j in range(dim)],
                               ]);

    def time_elapsed(self):
        return self.state[0][0];

    def objectsToState(self):
        n = len(objects);
        dim = len(objects[0].position);
        self.state = np.array([[t],
                               [[objects[i].position[j] for i in range(n)] for j in range(dim)],
                               [[objects[i].dirVec[j] for i in range(n)] for j in range(dim)],
                               ]);

    def stateToObjects(self):
        for i in range(len(objects)):
            objects[i].position = self.state[1][i];
            objects[i].dirVec = self.state[2][i];

    def velocityVector(self, i = None, name = None):
        if (i is None):
            if (name is not None):
                for j in range(len(self.objects)):
                    if (self.objects[j].name == name):
                        i = j;
            else:
                return np.zeroes(3);
            
        dist = [(fy.dist(self.objects[i], self.objects[j]) if (i != j) else 0.1) for j in range(len(objects))];      
        p0 = self.objects[i].position;

        vec = np.zeroes(3);
        for j in range(self.objects):
            if (i != j):
                r = dist[j];
                if (r == 0):
                    r = 0.1;
                pj = self.objects[j];
                vec += (pj.gm * (pj.position - p0)) / (r * r * r);
        return vec;

    def updateState(self, state):
        self.state = state;
        


                          
