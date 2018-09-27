import fy
from RungeKuttaFehlberg import RungeKuttaFehlberg54

import numpy as np
import time
import math
import matplotlib.pyplot as plot
import matplotlib.animation as animation

class System:
    def __init__(self,
                objects,
                 G = 6.67408 * 10**(-11), t = 0, estimate = None, stepsize = 0.25, tol = 0.0001):
        self.GravConst = G;
        self.objects = objects;
        self.objLen = len(objects);
        self.dim = len(objects[0].position);
        self.h = stepsize;
        self.tol = tol;
        self.state = np.array([np.array([t] * 2),
                               np.array([objects[i].position[0] for i in range(self.objLen)]),
                               np.array([objects[i].position[1] for i in range(self.objLen)]),
                               np.array([objects[i].dirVec[0] for i in range(self.objLen)]),
                               np.array([objects[i].dirVec[1] for i in range(self.objLen)])]);
        #self.state = np.array([np.array([t]),
        #                       np.array([np.array([objects[i].position[j] for i in range(self.objLen)]) for j in range(self.dim)]),
        #                       np.array([np.array([objects[i].dirVec[j] for i in range(self.objLen)]) for j in range(self.dim)]),
        #                       ]);
        self.estimate = estimate;

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
        p = np.array([self.state[1], self.state[2]]).T;
        v = np.array([self.state[3], self.state[4]]).T;
        #state1 = np.array(self.state[1]).T;
        #state2 = np.array(self.state[2]).T;
        self.objects[0].position = p[0];
        self.objects[1].position = p[1];
        self.objects[0].dirVec = v[0];
        self.objects[1].dirVec = v[1];
        #for i in range(self.objLen):
        #    self.objects[i].position = state1[i];
        #    self.objects[i].dirVec = state2[i];

    def velocityVector(self, x, i = None, name = None):
        if (i is None):
            if (name is not None):
                for j in range(self.objLen):
                    if (self.objects[j].name == name):
                        i = j;
            else:
                return np.zeros(3);
        p = np.asarray(x[1]);
        dist = self.distPos(p);
        p0 = p[:, i];

        vec = np.zeros(3);
        for j in range(self.objLen):
            if (i != j):
                r = dist[i][j];     # Can be simplified, there is only one i.
                if (r == 0):
                    r = 0.1;
                vec += (self.objects[j].gm * (p[:, j] - p0)) / (r * r * r);
        return vec;

    def distPos(self, x):
        return [[self.__tempSum(x, n, m)**(1.5) for m in range(self.objLen)] for n in range(self.objLen)];

    def __tempSum(self, x, n, m):
        sum0 = 0;
        for i in range(1, 3):    #self.dim
            temp = x[i][m] - x[i][n];
            sum0 += temp * temp;
        return sum0;

    def updateState(self, state):
        self.state = state;
        
    #def ydot(self, x):
    #    x = x;
    #    v = np.array([self.velocityVector(x, i) for i in range(self.objLen)]);
    #    return np.array([[1], x[1], v.T]);

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
        px = x[1];
        py = x[2];
        vx = x[3];
        vy = x[4];
        dist = self.distPos(x);
        #dist = [[((px[m] - px[n])**2 + (py[m] - py[n])**2)**(1.5) for m in range(self.objLen)]
        #        for n in range(self.objLen)];
        return np.array([np.ones(2), x[3], x[4], self.force(px, dist, Gm), self.force(py, dist, Gm)]);

    def rk_safestep(self):
        W, E = self.rk_step(self.h);
        # Check if the error is tolerable
        if (not (E < self.tol)):
            # Try to adjust the optimal step length
            if (E == 0):
                s = 2;
            else:
                s = math.pow(self.tol * self.h / (2 * E), 0.25);
            self.h = s * self.h;

            W, E = self.rk_step(self.h);
        # If the error is still not tolerable
        counter = 0;
        while (not (E < self.tol)):
            # Try if dividing the steplength with 2 helps.
            self.h = self.h/2;
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
        self.h = s * self.h;

        #Update self
        return W, E
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

    """def reshape_W(self, W):
        res = [W[0][0]];
        p_T = np.array(W[1]).T
        for p in p_T:
            for p_x in p:
                res.append(p_x);
        v_T = np.array(W[2]).T
        for v in v_T:
            for v_x in v:
                res.append(v_x);
        return np.array(res);

    def reshapeBack_W(self, W_out):
        l = self.objLen;
        d = self.dim;
        res = [[W_out[0]]];
        temp1 = [];
        w_x = 1;
        for j in range(l):
            temp2 = [];
            for i in range(w_x, w_x + d):
                temp2.append(W_out[i]);
                w_x += 1;
            temp1.append(temp2);
        res.append(np.array(temp1).T.tolist());
        temp1 = [];
        for j in range(l):
            temp2 = [];
            for i in range(w_x, w_x + d):
                temp2.append(W_out[i]);
                w_x += 1;
            temp1.append(temp2);
        res.append(np.array(temp1).T.tolist());
        return np.array(res);"""

    def step(self, rate = 3):
        W, E = self.rk_safe_step();
        self.updateState(W);
        for i in range(rate - 1):
            W, E = self.rk_safe_step();
            self.updateState(W);
        self.updateState(W);
        self.stateToObjects();
        return W, E;
                          
