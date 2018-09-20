import time
import matplotlib.animation as animation
import matplotlib.pyplot as plot
import numpy as np
from SaturnV import SaturnV
from RungeKuttaFehlberg import RungeKuttaFehlberg54


class Ascension(object):
    def __init__(self, init_state, grav_const, mass, planets):
        self.grav_const = grav_const
        self.mass = mass
        self.planets = planets
        self.state = np.array(init_state)

    def position(self):
        return [(self.state[1][n], self.state[2][n]) for n in range(self.planets)]

    def time_elapsed(self):
        return self.state[0][0]

    def step(self, h):
        print(saturn_v.get_mass(h))
        self.mass[1] = saturn_v.get_mass(h)

        x = self.state
        s1 = self.ydot(x)
        s2 = self.ydot(np.array([(np.array(f) * h) for f in s1]) + x)
        self.state = x + h * (s1 + s2) / 2

    def force(self, p, dist, gravity_mass):
        result = []
        for n in range(self.planets):
            temp_sum = 0
            for m in range(self.planets):

                if n != m and dist[n][m] > earth_radius:

                    if m == 2:
                        temp_sum += saturn_v.get_force(self.time_elapsed())

                    temp_sum += (gravity_mass[m] * (p[m] - p[n])) / (dist[n][m]**3)
                elif n != m and dist[n][m] < earth_radius:
                    # Rocket has crashed into the planet. Force and velocity is removed
                    self.state[3][1] = 0
                    self.state[4][1] = 0

            result.append(temp_sum)
        return np.array(result)

    def ydot(self, x):
        gravity_mass = self.mass * self.grav_const
        px = x[1]
        py = x[2]
        vx = x[3]
        vy = x[4]
        dist = [[((px[m] - px[n]) ** 2 + (py[m] - py[n]) ** 2) ** 0.5 for m in range(self.planets)] for n in range(self.planets)]

        return np.array([np.array([1]), vx, vy, self.force(px, dist, gravity_mass), self.force(py, dist, gravity_mass)])


earth_radius = 6.371 * 10**3
earth_mass = 5.972 * 10**24
rocket_mass = 2.97 * 10**6
grav_const = 6.67408 * 10**-11

# init_state is [t0,x0,y0,vx0,vy0]
planet = [0, 0, 0, 0, 0]
rocket = [0, 0, earth_radius + earth_radius/10, 0, 150000]
init = np.array([
    np.array([0.0]),
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
axes = fig.add_subplot(111, aspect="equal", autoscale_on=False, xlim=(-5*earth_radius, 5*earth_radius), ylim=(-5*earth_radius, 5*earth_radius))

line1, = axes.plot([], [], "o-g", lw=2, ms=earth_radius/125)  # A blue planet
line2, = axes.plot([], [], "2-r", lw=2, ms=earth_radius/1000)  # A red rocket ship
time_text = axes.text(0.02, 0.95, "", transform=axes.transAxes)


def init():
    # initialize animation
    line1.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    return line1, line2, time_text


def animate(i):
    # perform animation step
    global orbit, dt
    for i in range(1):
        orbit.step(dt)
    pos = orbit.position()
    line1.set_data(*pos[0])
    line2.set_data(*pos[1])

    # Scale earth size
    # left, right = plot.xlim()
    # line1.set_markersize(0.0008*(right-left))
    # print(right-left)

    time_text.set_text('time = %.1f' % orbit.time_elapsed())
    return line1, line2, time_text


# choose the interval based on dt and the time to animate one step
# Take the time for one call of the animate.
t0 = time.time()
animate(0)
t1 = time.time()

delay = 1000 * dt - (t1 - t0)

anim = animation.FuncAnimation(fig,             # figure to plot in
                               animate,         # function that is called on each frame
                               frames=3000,        # total number of frames
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
