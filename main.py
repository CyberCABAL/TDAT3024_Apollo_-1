from fy import CelestialObject
from orbitalsys import System

import numpy as np
import time
import fy
import matplotlib.pyplot as plot
import matplotlib.animation as animation

dt = 1./30 # 30 frames per second
tol = 3e-14;
sys = System([CelestialObject([0., 0.], [0., 0.], 5.9722 * 10**24, 6378100, "Terra"),
              CelestialObject([362600000., 0.], [0., 1078.2], 7.34767309 * 10**22, 1737000, "Luna")], stepsize = dt, tol = tol);
winDimention = 362600000. * 2;
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

def main():
    # The figure is set
    fig = plot.figure();
    axes = fig.add_subplot(111, aspect="equal", autoscale_on=False, xlim=(-winDimention, winDimention), ylim=(-winDimention, winDimention));
    line1, = axes.plot([], [], "o-b", lw=2); # Terra
    line2, = axes.plot([], [], "o-k", lw=2); # Luna
    time_text = axes.text(0.02, 0.95, "", transform=axes.transAxes);
    dist_text = axes.text(0.02, 0.90, "", transform=axes.transAxes);

    def init():
        #initialize animation
        line1.set_data([], []);
        line2.set_data([], []);
        time_text.set_text('');
        dist_text.set_text("");
        return line1, line2, time_text, dist_text;

    def animate(i):
        #perform animation step
        global sys;
        sys.step(4);
        line1.set_data(*sys.objects[0].position);
        line2.set_data(*sys.objects[1].position);
        time_text.set_text('time = %.1f' % (sys.time_elapsed()/60/60/24) + " days");
        dist_text.set_text("distance = %.1f" % fy.dist(sys.objects[0], sys.objects[1]));
        return line1, line2, time_text, dist_text;

    # choose the interval based on dt and the time to animate one step
    # Take the time for one call of the animate.
    t0 = time.time();
    animate(0);
    t1 = time.time();

    delay = 1000 * dt - (t1 - t0);

    anim=animation.FuncAnimation(fig,        # figure to plot in
                            animate,    # function that is called on each frame
                            frames=3000, # total number of frames 
                            interval=delay, # time to wait between each frame.
                            repeat=False,
                            blit=True, 
                            init_func=init # initialization
                            );

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    #anim.save('orbit.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

    plot.show();
    
if __name__ == "__main__":
    main();
