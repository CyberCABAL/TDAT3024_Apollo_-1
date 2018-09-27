from numpy import sqrt
import time
from RungeKuttaFehlberg import RungeKuttaFehlberg54
import numpy as np
import math

import matplotlib.pyplot as plot
import matplotlib.animation as animation

class Orbit:
    def __init__(self,
                 init_state = [[0, 0, 0], [0, 1, 2], [0, 1, 2], [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]],
                 G=1,
                 mass = [1, 1, 1],
                 planets = 3,
                 stepsize = 0.15,
                 tolerance = 05e-15):
        self.GravConst = G;
        self.m = mass;
        self.planets = planets;
        self.state = np.array(init_state);
        self.h = stepsize
        self.tol = tolerance
    
    def position(self):
        #compute the current x, y positions of the pendulum arms
        return [((self.state[1][n], self.state[2][n])) for n in range(self.planets)];

    def time_elapsed(self):
        return self.state[0][0];

    def rk_safestep(self):
        w, e = self.rk_step(self.h);
        # Check if the error is tolerable
        if (not (e < self.tol)):
            # Try to adjust the optimal step length
            if (e == 0):
                s = 2;
            else:
                s = math.pow(self.tol * self.h / (2 * e), 0.25);
            self.h = s * self.h;

            w, e = self.rk_step(self.h);
        # If the error is still not tolerable
        counter = 0;
        while (not (e < self.tol)):
            # Try if dividing the steplength with 2 helps.
            self.h = self.h/2;
            w, e = self.rk_step(self.h);
            counter += 1;
            if (counter > 10):
                print("System is unreliable, terminating.")
                sys.exit(-1);

        #Adjust step
        if (e == 0):
            s = 2;
        else:
            s = math.pow(self.tol * self.h / (2 * e), 0.25);
        self.h = s * self.h;

        #Update self
        self.state = w

    def rk_step(self, h):
        x = self.state
        s1 = self.ydot(x)
        s2 = self.ydot(x+h*1/4*s1)
        s3 = self.ydot(x+h*(3/32*s1 + 9/32*s2))
        s4 = self.ydot(x+h*(1932/2197*s1 + -7200/2197*s2 + 7296/2197*s3))
        s5 = self.ydot(x+h*(439/216*s1 + -8*s2 + 3680/513*s3 + -845/4104*s4))
        s6 = self.ydot(x+h*(-8/27*s1 + 2*s2 + -3544/2565*s3 + 1859/4104*s4 + -11/40*s5))

        w = x + h* (25/216*s1 + 1408/2565*s3 + 2197/4104*s4 + -1/5*s5)
        z = x + h* (16/135*s1 + 6656/12825*s3 + 28561/56430*s4 + -9/50*s5 + 2/55*s6)

        e = np.linalg.norm(w - z, 2) / np.linalg.norm(w, 2)

        print(w)
        self.state = w

        return w, e

    def step(self, h):
        #Uses the trapes method to calculate the new state after h seconds.
        x = self.state;
        s1 = self.ydot(x);
        s2 = self.ydot(np.array([(np.array(f) * h) for f in s1]) + x)
        self.state = x + h * (s1 + s2) / 2;
        print(self.state)
        
    def eulerstep(self, h):
        #Uses the euler method to calculate the new state after h seconds.
        x = self.state;
        self.state = x + h * self.ydot(x);

    def force(self, p, dist, Gm):
        result = [];
        for n in range(self.planets):
            tempSum = 0;
            for m in range(self.planets):
                if (n != m and dist[n][m] != 0):
                    tempSum += (Gm[m] * (p[m] - p[n])) / (dist[n][m]);
            result.append(tempSum);
        return np.array(result);
    
    def ydot(self, x):
        '''''
        y1 =
        '''
        Gm = self.m * self.GravConst;
        px = x[1];
        py = x[2];
        vx = x[3];
        vy = x[4];
        dist = [[((px[m] - px[n])**2 + (py[m] - py[n])**2)**(1.5) for m in range(self.planets)]
                for n in range(self.planets)];
                
        return np.array([np.ones(self.planets), vx, vy, self.force(px, dist, Gm), self.force(py, dist, Gm)]);

# make an Orbit instance
orbit = Orbit(init_state=init, G=G, mass=m, planets=2);
dt = 1./30 # 30 frames per second

# The figure is set
fig = plot.figure();
axes = fig.add_subplot(111, aspect="equal", autoscale_on=False, xlim=(-1000, 1000), ylim=(-1000, 1000))

line1, = axes.plot([], [], "o-b", ms=12.7); # A green planet
line2, = axes.plot([], [], marker='o', color='tab:gray', ms=1.7); # A red planet
# line3, = axes.plot([], [], "o-g", lw=2); # A blue planet
time_text = axes.text(0.02, 0.95, "", transform=axes.transAxes);

def init():
    #initialize animation
    line1.set_data([], []);
    line2.set_data([], []);
    # line3.set_data([], []);
    time_text.set_text('');
    return line1, line2, time_text;

def animate(i):
    #perform animation step
    global orbit, dt;
    for i in range(1):
        orbit.rk_safestep();
    pos = orbit.position();
    line1.set_data(*pos[0]);
    line2.set_data(*pos[1]);
    # line3.set_data(*pos[2]);
    time_text.set_text('time = %.1f' % orbit.time_elapsed());
    return line1, line2, time_text;

# choose the interval based on dt and the time to animate one step
# Take the time for one call of the animate.
t0 = time.time();
animate(0);
t1 = time.time();

delay = 1000 * dt - (t1 - t0);

anim=animation.FuncAnimation(fig,        # figure to plot in
                        animate,    # function that is called on each frame
                        frames = 3000, # total number of frames 
                        interval = delay, # time to wait between each frame.
                        repeat = False,
                        blit = True, 
                        init_func=init # initialization
                        );

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('orbit.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plot.show();
