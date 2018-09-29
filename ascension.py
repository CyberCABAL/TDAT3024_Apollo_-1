import time
import matplotlib.animation as animation
import matplotlib.pyplot as plot
import numpy as np
import fy
import math
from fy import CelestialObject
from SaturnV import SaturnV


class Ascension(object):
    def __init__(self, init_state, grav_const, mass, planets,
                 stepsize=0.15,
                 tolerance=05e-14):
        self.grav_const = grav_const
        self.mass = mass
        self.planets = planets
        self.state = np.array(init_state)
        self.h = stepsize
        self.tol = tolerance

    def position(self):
        return [(self.state[1][n], self.state[2][n]) for n in range(self.planets)]

    def time_elapsed(self):
        return self.state[0][0]

    def rk_safestep(self):
        w, e = self.rk_step(self.h);
        # Check if the error is tolerable
        if (not (e < self.tol)):
            # Try to adjust the optimal step length
            if (e == 0):
                s = 2;
            else:
                s = math.pow(self.tol * self.h / (2 * e), 0.25);
                if(s>2):
                    s=2
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
            if (s > 2):
                s = 2
        self.h = s * self.h;

        #Update self
        self.state = w

    def rk_step(self, h):
        x = self.state
        #print('\n\nsteplength: '+str(h))
        #print('s1')
        s1 = self.ydot(x)
        #print('s2')
        s2 = self.ydot(x+h*1/4*s1)
        #print('s3')
        s3 = self.ydot(x+h*(3/32*s1 + 9/32*s2))
        #print('s4')
        s4 = self.ydot(x+h*(1932/2197*s1 + -7200/2197*s2 + 7296/2197*s3))
        #print('s5')
        s5 = self.ydot(x+h*(439/216*s1 + -8*s2 + 3680/513*s3 + -845/4104*s4))
        #print('s6')
        s6 = self.ydot(x+h*(-8/27*s1 + 2*s2 + -3544/2565*s3 + 1859/4104*s4 + -11/40*s5))

        w = x + h* (25/216*s1 + 1408/2565*s3 + 2197/4104*s4 + -1/5*s5)
        z = x + h* (16/135*s1 + 6656/12825*s3 + 28561/56430*s4 + -9/50*s5 + 2/55*s6)

        e = np.linalg.norm(w - z, 2) / np.linalg.norm(w, 2)

        # print(w)
        # self.state = w

        return w, e

    def step(self, h):
        self.mass[1] = saturn_v.get_mass(self.time_elapsed())

        x = self.state
        s1 = self.ydot(x)
        s2 = self.ydot(np.array([(np.array(f) * h) for f in s1]) + x)
        self.state = x + h * (s1 + s2) / 2

    #Actually acceleration
    def force(self, p, dist):
        G = self.grav_const
        mass = self.mass
        result = []
        for n in range(self.planets):
            temp_sum = 0

            #Rocket motor
            if n == 1:
                F = saturn_v.get_force(self.time_elapsed())
                temp_sum += F / saturn_v.get_mass(self.time_elapsed())

                #temp_sum += G * (mass[0] * (p[0] - p[n])) / (dist[n][0] ** 3)
                #print('f acc: '+str(temp_sum) + '\n')

            #Gravitational forces
            
            for m in range(self.planets):
                if n != m and dist[n][m] > earth_radius:
                    F = G*(mass[m]*(p[m] - p[n])) / (dist[n][m]**3)
                    temp_sum += F
                elif n != m and dist[n][m] < earth_radius:
                    # Rocket has crashed into the planet. Velocity is removed from rocket
                    self.state[3][1] = 0
                    self.state[4][1] = 0

                    # F = G * (mass[m] * mass[n] * (p[m] - p[n])) / (dist[n][m] ** 3)
                    # temp_sum += -F / mass[n]

            result.append(temp_sum)
        return np.array(result)

    def ydot(self, x):
        px = x[1]
        py = x[2]
        vx = x[3]
        vy = x[4]
        dist = [[((px[m] - px[n]) ** 2 + (py[m] - py[n]) ** 2) ** 0.5 for m in range(self.planets)] for n in range(self.planets)]

        return np.array([np.ones(self.planets), vx, vy, self.force(px, dist), self.force(py, dist)])


earth_radius = 6.371 * 10**6
earth_mass = 5.972 * 10**24
rocket_mass = 2.97 * 10**6
grav_const = 6.67408 * 10**-11

# init_state is [t0,x0,y0,vx0,vy0]
planet = [0, 0, 0, 0, 0]
rocket = [0, 0, earth_radius+10, 0, 0]
init = np.array([
    np.array([0.0, 0]),
    np.array([planet[1], rocket[1]]),
    np.array([planet[2], rocket[2]]),
    np.array([planet[3], rocket[3]]),
    np.array([planet[4], rocket[4]])
])
# G = 6.67408 * 10**(-11);	# Real G
mass = np.array([earth_mass, rocket_mass])
planets = 2

# make an Orbit instance
orbit = Ascension(init, grav_const, mass, planets)
dt = 1/300  # 300 frames per second

saturn_v = SaturnV()

# The figure is set
fig = plot.figure()
viewsize = 5
axes = fig.add_subplot(111, aspect="equal", autoscale_on=False, xlim=(-viewsize*earth_radius, viewsize*earth_radius), ylim=(-viewsize*earth_radius, viewsize*earth_radius))

trail, = axes.plot([], [], "o-b", lw=1, ms=earth_radius/100000000, label='trail')  # A blue dotted trail
line1, = axes.plot([], [], "o-g", lw=0, ms=earth_radius/125000, label='Earth')  # A blue planet
line2, = axes.plot([], [], "2-r", lw=0, ms=earth_radius/1000000, label='Rocket')  # A red rocket ship
time_text = axes.text(0.02, 0.95, "", transform=axes.transAxes)
posx_text = axes.text(0.02, 0.90, "", transform=axes.transAxes)
posy_text = axes.text(0.02, 0.85, "", transform=axes.transAxes)
legend = axes.legend(loc='lower right')
for legend_handle in legend.legendHandles:
    legend_handle._legmarker.set_markersize(6)

trailx = [rocket[1]]
traily = [rocket[2]]
framecounter = 1;

def init():
    # initialize animation
    line1.set_data([], [])
    line2.set_data([], [])
    trail.set_data([], [])
    time_text.set_text('')
    posx_text.set_text('')
    posy_text.set_text('')
    return line1, line2, time_text, posx_text, posy_text, trail


def animate(i):
    # perform animation step
    global orbit, dt, framecounter

    # while(orbit.time_elapsed() < framecounter/30):
    for i in range(1):
        orbit.rk_safestep()
    # framecounter+= 1

    # for i in range(10):
    #     pass
    pos = orbit.position()
    trailx.append(pos[1][0])
    traily.append(pos[1][1])
    line1.set_data(*pos[0])
    line2.set_data(*pos[1])
    trail.set_data(trailx,traily)

    # Scale earth size
    left, right = plot.xlim()
    line1.set_markersize(earth_radius/((right-left)/510))
    # print(right-left)

    time_text.set_text('time = %.2f' % orbit.time_elapsed())
    posx_text.set_text('posx =  %.3e' % orbit.position()[1][0])
    posy_text.set_text('posy =  %.3e' % orbit.position()[1][1])
    return line1, line2, time_text, posx_text, posy_text, trail


# choose the interval based on dt and the time to animate one step
# Take the time for one call of the animate.
t0 = time.time()
animate(0)
t1 = time.time()

delay = 1000 * dt - (t1 - t0)

anim = animation.FuncAnimation(fig,             # figure to plot in
                               animate,         # function that is called on each frame
                               frames=30000,        # total number of frames
                               interval=delay,  # time to wait between each frame.
                               repeat=False,
                               blit=False,
                               init_func=init   # initialization
                               )

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
# anim.save('orbit.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plot.show()
