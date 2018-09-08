
from numpy import sqrt
import time

import numpy as np
import scipy.integrate as integrate

import matplotlib.pyplot as plot
import matplotlib.animation as animation

#G = 6.67408 * 10**(-11);
G = 1.0;
m = np.array([1, 1, 1]);
#initXY = [[2, 2], [0, 0], [-2, -2]];
#initdXdY = [[0.2, -0.2], [0, 0], [-0.2, 0.2]];
init = np.array([
    np.array([0.0]),
    np.array([-0.970, 0.970, 0.0]),
    np.array([0.243, -0.243, 0.0]),
    np.array([-0.466, -0.466, 0.932]),
    np.array([-0.433, -0.433, 0.866])]);

class Orbit:
    """
    
    Orbit Class

    init_state is [t0,x0,vx0,y0,vx0],
    where (x0,y0) is the initial position
    , (vx0,vy0) is the initial velocity
    and t0 is the initial time
    """
    def __init__(self,
                 init_state = [[0], [0, 1, 2], [0, 1, 2], [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]],
                 G=1,
                 mass = [1, 1, 1],
                 planets = 3):
        self.GravConst = G
        self.m = mass;
        self.planets = planets;
        self.state = np.array(init_state);
    
    def position(self):
        #compute the current x, y positions of the pendulum arms
        return [((self.state[1][n], self.state[2][n])) for n in range(self.planets)];

    def time_elapsed(self):
        return self.state[0][0]

    def step(self, h):
        """Uses the trapes method to calculate the new state after h seconds."""
        x = self.state
        s1 = self.ydot(x)
        s2 = self.ydot(np.array([(np.array(f) * h) for f in s1]) + x)
        self.state = x + h * (s1 + s2) / 2
        #self.state = [((np.array(s1[n] + s2[n]) / 2)) for n in range(len(s1))] + x;
        
    def eulerstep(self, h):
        """Uses the euler method to calculate the new state after h seconds."""
        x = self.state;
        self.state = x + h * self.ydot(x);
        #dot = self.ydot(x);
        #self.state = np.array([(np.array(f) * h) for f in dot]) + x;

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
        Gm = self.m * self.GravConst;
        px = x[1];
        py = x[2];
        vx = x[3];
        vy = x[4];
        dist = [[((px[m] - px[n])**2 + (py[m] - py[n])**2)**(1.5) for m in range(self.planets)]
                for n in range(self.planets)];
        #print(dist);
                
        return np.array([np.array([1]), vx, vy, self.force(px, dist, Gm), self.force(py, dist, Gm)]);

# make an Orbit instance
orbit = Orbit(init, G, m, 3);
dt = 1./30 # 30 frames per second

# The figure is set
fig = plot.figure()
axes = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-10, 10), ylim=(-10, 10))

line1, = axes.plot([], [], 'o-g', lw=2) # A green planet
line2, = axes.plot([], [], 'o-b', lw=2) # A blue planet
line3, = axes.plot([], [], 'o-r', lw=2) # A red planet
time_text = axes.text(0.02, 0.95, '', transform=axes.transAxes)

def init():
    """initialize animation"""
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    time_text.set_text('')
    return line1, line2, line3, time_text

def animate(i):
    """perform animation step"""
    global orbit, dt
    for i in range(10):
        orbit.step(dt)
    pos = orbit.position();
    line1.set_data(*pos[0])
    line2.set_data(*pos[1])
    line3.set_data(*pos[2])
    time_text.set_text('time = %.1f' % orbit.time_elapsed())
    return line1, line2, line3, time_text

# choose the interval based on dt and the time to animate one step
# Take the time for one call of the animate.
t0 = time.time()
animate(0)
t1 = time.time()

delay = 1000 * dt - (t1 - t0)

anim=animation.FuncAnimation(fig,        # figure to plot in
                        animate,    # function that is called on each frame
                        frames=3000, # total number of frames 
                        interval=delay, # time to wait between each frame.
                        repeat=False,
                        blit=True, 
                        init_func=init # initialization
                        )

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('orbit.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plot.show()

"""
K = 1: Planetene beveger seg ut av mønsteret etter hvert, og systemet blir ødelagt. Etter hvert kommer planetene så nær hverandre at de blir kastet ut.
K = 2: Det skjer liten endring over tid, men det virker som om hele systemet flytter seg litt i rommet, så regner med at det vil gå galt etter lang tid.
K = 3 - 5: Lite endring, ganske stabil. Resultat ble mer eller mindre nøyaktig likt.
"""
