from fy import CelestialObject
from fy import Rocket
from orbitalsys import System

import numpy as np
import time
import fy
import matplotlib.pyplot as plot
import matplotlib.animation as animation

dt = 1./30 # 30 frames per second
tol = 3e-14
x_0 = -6378100 * 2/3
dy_0 = -13.098
luna_distance = 362000000.
terra_r = 6378100
luna_r = 1737000
sys = System([CelestialObject([x_0, 0.], [0., dy_0], 5.9722 * 10**24, terra_r, "Terra"),
              CelestialObject([luna_distance + x_0, 0.], [0., 1078.2 + dy_0], 7.34767309 * 10**22, luna_r, "Luna"),
              Rocket([x_0, terra_r + 10], [0.05, 1], 2.97 * 10**6, 1, "Saturn V")
              ], stepsize=dt, tol=tol)
winDimention = luna_distance * 3/2

t_colour = [0.1, 0.33, 0.8]
l_colour = [0.5, 0.47, 0.55]
bk_colour = [0.17, 0.2, 0.3]
tr_colour = [0.88, 0.6, 0.15]

"""
p'1 = v1
v'1 = Gm2(p2 - p1)/r**3_12
p'2 = v2
v'2 = Gm1(p1 - p2)/r**3_12

p1(0) = [0,0,0]
p2(0) = [362600., 0., 0.]
p'1(0) = [0., 0., 0.]
p'2(0) = [0., 1078.2, 0.]
"""
line1 = None
line2 = None
line3 = None
trail = None

posx_text = None
posy_text = None

trailx = [sys.objects[2].position[0]]
traily = [sys.objects[2].position[1]]

def main():
    # The figure is set
    fig = plot.figure()
    axes = fig.add_subplot(111, aspect="equal", autoscale_on=False, xlim=(-winDimention, winDimention), ylim=(-winDimention, winDimention))
    axes.set_facecolor(bk_colour)
    # line1, = axes.plot([], [], "o-b", lw=2) # Terra
    # line2, = axes.plot([], [], "o-k", lw=2) # Luna
    global line1, line2, line3, trail
    line1 = plot.Circle((0, 0), terra_r, color=t_colour)
    line2 = plot.Circle((luna_distance, 0), luna_r, color=l_colour)
    line3, = axes.plot([], [], "d-w", lw=0, ms=terra_r/5000000, label="Rocket")
    trail, = axes.plot([], [], color=tr_colour, marker="o", lw=0.5, ms=terra_r/10000000000, label='Trail', alpha=0.85, ls="--")
    axes.add_artist(line1)
    axes.add_artist(line2)
    posx_text = axes.text(0.02, 0.80, "", transform=axes.transAxes, color="w")
    posy_text = axes.text(0.02, 0.85, "", transform=axes.transAxes, color="w")
    time_text = axes.text(0.02, 0.95, "", transform=axes.transAxes, color="w")
    dist_text = axes.text(0.02, 0.90, "", transform=axes.transAxes, color="w")

    def init():
        #initialize animation
        #line1
        #line2
        line3.set_data([], [])
        time_text.set_text("")
        dist_text.set_text("")
        posx_text.set_text("")
        posy_text.set_text("")
        trail.set_data([], [])
        return line1, line2, line3, time_text, dist_text, trail, posx_text, posy_text

    def animate(i):
        #perform animation step
        global sys, line1, line2, line3, trail
        sys.step(6)    #originally 4
        #line1.set_data(*sys.objects[0].position)
        #line2.set_data(*sys.objects[1].position)
        line1.remove()
        line2.remove()
        del line1
        del line2
        line1 = plot.Circle((sys.objects[0].position[0], sys.objects[0].position[1]), 6378100, color=t_colour)
        line2 = plot.Circle((sys.objects[1].position[0], sys.objects[1].position[1]), 1737000, color=l_colour)
        axes.add_artist(line1)
        axes.add_artist(line2)
        trailx.append(sys.objects[2].position[0])
        traily.append(sys.objects[2].position[1])
        trail.set_data(trailx, traily)
        line3.set_data((sys.objects[2].position[0], sys.objects[2].position[1]))
        #line1.xy = sys.objects[0].position
        #line2.xy = sys.objects[1].position
        time_text.set_text('time = %.1f' % (sys.time_elapsed()/86400) + " days")
        dist_text.set_text("distance = %.1f" % fy.dist(sys.objects[0], sys.objects[1]))
        posx_text.set_text('posx =  %.3e' % sys.objects[2].position[0])
        posy_text.set_text('posy =  %.3e' % sys.objects[2].position[1])
        return line1, line2, line3, time_text, dist_text, trail, posx_text, posy_text

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
                            blit=False,
                            init_func=init # initialization
                            )

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    #anim.save('orbit.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

    plot.show()
    
if __name__ == "__main__":
    main()
