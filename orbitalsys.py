import fy

import numpy as np
import time
import math
import matplotlib.pyplot as plot
import matplotlib.animation as animation

class System:
    
    def __init__(self,
                objects,
                 G = 6.67408 * 10**(-11), t = 0, stepsize = 0.25, tol = 0.0001):
        self.GravConst = G;
        self.objects = objects;
        self.objLen = len(objects);
        self.dim = 2;#len(objects[0].position);
        self.h = stepsize;
        self.tol = tol;

        print("dim:", self.dim);
        
        st = [np.array([t] * self.dim)];
        for j in range(self.dim):
            st.append(np.array([objects[i].position[j] for i in range(self.objLen)]));
        for j in range(self.dim):
            st.append(np.array([objects[i].dirVec[j] for i in range(self.objLen)]));
        self.state = np.array(st);

    def time_elapsed(self):
        return self.state[0][0];

    def objectsToState(self):
        self.state = np.array([self.state[0],
                               [[self.objects[i].position[j] for i in range(self.objLen)] for j in range(self.dim)],
                               [[self.objects[i].dirVec[j] for i in range(self.objLen)] for j in range(self.dim)],
                               ]);

    def addObject(self, obj):
        self.objects.append(obj);
        self.objLen += 1;

    def stateToObjects(self):
        p = [];
        for i in range(1, 1 + self.dim):
            p.append(self.state[i]);
        p = np.array(p).T;
        v = [];
        for i in range(-self.dim, 0):
            v.append(self.state[i]);
        v = np.array(v).T;
        for i in range(self.dim):
            self.objects[i].position = p[i];
            self.objects[i].dirVec = v[i];

    def distPos(self, x):
        return [[self.__tempSum(x, n, m)**(1.5) for m in range(self.objLen)] for n in range(self.objLen)];

    def __tempSum(self, x, n, m):
        sum0 = 0;
        for i in range(1, self.dim + 1):
            temp = x[i][m] - x[i][n];
            sum0 += temp * temp;
        return sum0;

    def updateState(self, state):
        self.state = state;

    def force(self, p, dist, Gm):
        result = [];
        for n in range(self.objLen):
            tempSum = 0;
            for m in range(self.objLen):
                if (n != m and dist[n][m] != 0):
                    tempSum += (Gm[m] * (p[m] - p[n])) / (dist[n][m]);
            result.append(tempSum);
        return np.array(result);

    def ydot(self, x):
        Gm = np.array([o.mass for o in self.objects]) * self.GravConst;
        dist = self.distPos(x);
        res = [np.ones(self.dim)];
        for i in range(-self.dim, 0):
            res.append(x[i]);
        for i in range(1, 1 + self.dim):
            res.append(self.force(x[i], dist, Gm));
        return np.array(res);
        #return np.array([np.ones(self.dim), x[3], x[4], self.force(x[1], dist, Gm), self.force(x[2], dist, Gm)]);

    def rk_safestep(self):
        W, E = self.rk_step(self.h);
        # Check if the error is tolerable
        if (not (E < self.tol)):
            # Try to adjust the optimal step length
            if (E == 0):
                s = 2;
            else:
                s = math.pow(self.tol * self.h / (2 * E), 0.25);
            self.h *= s;

            W, E = self.rk_step(self.h);
        # If the error is still not tolerable
        counter = 0;
        while (not (E < self.tol)):
            # Try if dividing the steplength with 2 helps.
            self.h /= 2;
            W, E = self.rk_step(self.h);
            counter += 1;
            if (counter > 10):
                print("System is unreliable, terminating.")
                sys.exit(-1);

        #Adjust step
        if (E == 0):
            s = 2;
        else:
            s = math.pow(self.tol * self.h / (2 * E), 0.25);
        self.h *= s;

        return W, E;
        #self.stateToObjects();

    def rk_step(self, h):
        x = self.state
        
        s1 = self.ydot(x)
        s2 = self.ydot(x+h*0.25*s1)
        s3 = self.ydot(x+h*(3/32*s1 + 9/32*s2))
        s4 = self.ydot(x+h*(1932/2197*s1 + -7200/2197*s2 + 7296/2197*s3))
        s5 = self.ydot(x+h*(439/216*s1 + -8*s2 + 3680/513*s3 + -845/4104*s4))
        s6 = self.ydot(x+h*(-8/27*s1 + 2*s2 + -3544/2565*s3 + 1859/4104*s4 + -11/40*s5))

        w = x + h* (25/216*s1 + 1408/2565*s3 + 2197/4104*s4 + -1/5*s5)
        z = x + h* (16/135*s1 + 6656/12825*s3 + 28561/56430*s4 + -9/50*s5 + 2/55*s6)

        e = np.linalg.norm(w - z, 2) / np.linalg.norm(w, 2)

        return w, e

    def step(self, rate = 3):
        W, E = self.rk_safestep();
        self.updateState(W);
        for i in range(rate - 1):
            W, E = self.rk_safestep();
            self.updateState(W);
        self.stateToObjects();
        return W, E;
                          
