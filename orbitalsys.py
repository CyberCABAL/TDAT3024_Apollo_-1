import fy

import numpy as np
import time
import math
import matplotlib.pyplot as plot
import matplotlib.animation as animation

crash_index = None;

class System:
    
    def __init__(self,
                objects,
                G = 6.67408 * 10**(-11), t = 0, stepsize = 0.25, tol = 0.0001, r_index = None):
        self.GravConst = G;
        self.objects = objects;
        self.objLen = len(objects);
        self.dim = len(objects[0].position);    #2
        self.time_dim = 1    #max(len(objects[0].position), len(objects));
        self.h = stepsize;
        self.tol = tol;
        self.r_index = r_index;

        #print("dim:", self.dim);
        
        st = [np.array([t] * self.time_dim)];
        # [[time], [x0, x1, ..., x_n], [y0, y1, ..., y_n], ..., [dx0, dx1, ..., dx_n], [dy0, dy1, ..., dy_n], ...]
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
        for i in range(self.time_dim, self.time_dim + self.dim):
            p.append(self.state[i]);
        p = np.array(p).T;
        v = [];
        for i in range(-self.dim, 0):
            v.append(self.state[i]);
        v = np.array(v).T;
        for i in range(self.objLen):
            self.objects[i].position = p[i];
            self.objects[i].dirVec = v[i];

    def distPos(self, x):
        return [[self.__tempSum(x, n, m)**(0.5) for m in range(self.objLen)] for n in range(self.objLen)];

    def __tempSum(self, x, n, m):
        sum0 = 0;
        for i in range(self.time_dim, self.dim + self.time_dim):
            temp = x[i][m] - x[i][n];
            sum0 += temp * temp;
        return sum0;

    def updateState(self, state):
        self.state = state;

    def a_G(self, p, dist, dist3, Gm, r_index):
        result = [];
        for n in range(self.objLen):
            tempSum = 0;
            if (not self.objects[n].stop):
                for m in range(self.objLen):
                    if (not self.objects[m].stop):
                        if (n != m and dist[n][m] > self.objects[m].r):
                            tempSum += (Gm[m] * (p[m] - p[n])) / (dist3[n][m]);
                        elif n != m and dist[n][m] < self.objects[m].r and n == r_index:
                            # Object n has crashed into object m. Velocity is removed.
                            self.objects[r_index].stop = True;
                            print("Crashed into", self.objects[m].name, "!");
                            global crash_index;
                            crash_index = m;
            result.append(tempSum);
        return np.array(result);

    def Σa(self, x, dist, r_index):
        Gm = np.array([o.mass for o in self.objects]) * self.GravConst;
        dist3 = [[dist[n][m] ** 3 for m in range(self.objLen)] for n in range(self.objLen)];

        if (r_index):
            rocket = self.objects[r_index];
        a = [self.a_G(x[i], dist, dist3, Gm, r_index) for i in range(self.time_dim, self.time_dim + self.dim)];
        if (r_index):
            if (not rocket.stop):
                a_R = rocket.a_R();
                a_A = rocket.a_Atmos(self.objects[0]);
                for i in range(self.dim):
                    a[i][r_index] += a_R[i] + a_A[i];
        return a;

    def ydot(self, x):
        index = self.r_index;
        if (index):
            self.objects[index].update(self.time_elapsed());
        #Gm = np.array([o.mass for o in self.objects]) * self.GravConst;
        dist = self.distPos(x);
        #print(dist)
        dist3 = [[dist[n][m] * dist[n][m] * dist[n][m] for m in range(self.objLen)] for n in range(self.objLen)];
        #print(dist3)
        res = [np.ones(self.time_dim)];
        for i in range(-self.dim, 0):
            res.append(x[i]);
        Σ = self.Σa(x, dist, index);
        for i in range(self.dim):
            res.append(Σ[i]);

        if (index):
            if (self.objects[index].stop):
                for i in range(-self.dim, 0):
                    res[i][index] = 0;
                for i in range(self.time_dim, self.time_dim + self.dim):
                    res[i][index] = res[i][crash_index];
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

        #w[0] = np.repeat(w[0], self.time_multiply);
        #z[0] = np.repeat(z[0], self.time_multiply);
        #print(w, z)
        e = np.linalg.norm(w - z, 2) / np.linalg.norm(w, 2)

        #w[0] = np.array([w[0][0]]);
        #z[0] = np.array([z[0][0]]);
        #print(w, z)

        return w, sum(e)/len(e);    #Uncertain about this.

    def step(self, rate = 3):
        for i in range(rate):
            W, E = self.rk_safestep();
            self.updateState(W);
            self.stateToObjects();
        if (self.r_index):
            if (self.objects[self.r_index].stop):
                self.objects[self.r_index].update(self.time_elapsed());
        return W, E;
                          
