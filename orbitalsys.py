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
        self.shape = np.array(objects).shape;
        self.state = np.array([[t],
                               [[objects[i].position[j] for i in range(self.objLen)] for j in range(self.dim)],
                               [[objects[i].dirVec[j] for i in range(self.objLen)] for j in range(self.dim)],
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
        state1 = np.array(self.state[1]).T;
        state2 = np.array(self.state[2]).T;
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
        for i in range(0, self.dim):
            temp = x[i][m] - x[i][n];
            sum0 += temp * temp;
        return sum0;

    def updateState(self, state):
        self.state = state;
        
    def ydot(self, x):
        x = self.reshapeBack_W(x);
        v = np.array([self.velocityVector(x, i) for i in range(self.objLen)]);
        return self.reshape_W(np.array([[1], x[1], v.T]));

    def reshape_W(self, W):
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
        return np.array(res);

    def step(self, rate = 3):
        W, E = self.estimate.safe_step(self.reshape_W(self.state));
        for i in range(rate - 1):
            W, E = self.estimate.safe_step(W);
        self.updateState(self.reshapeBack_W(W));
        self.stateToObjects();
        print(self.state);
        return W, E;

    #def step2(self, rate):
    #    r = None;
    #    for j in range(objLen):
    #       w = np.array([self.state[0], self.state[1][j], self.state[2][j]]);
    #        for i in range(rate):
    #            r = self.estimate[j].safeStep(w);
    #    return r;
                          
