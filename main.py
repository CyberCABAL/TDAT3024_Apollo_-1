from fy import CelestialObject
from orbitalsys import System
from RungeKuttaFehlberg import RungeKuttaFehlberg54

import numpy as np
import time
import matplotlib.pyplot as plot
import matplotlib.animation as animation
import matplotlib.patches as patch

sys = System([CelestialObject([0., 0., 0.], [0., 0., 0.], 59722. * 10**20, 6378100, "Terra"),
              CelestialObject([362600., 0., 0.], [0., 1078.2, 0.], 734767309. * 10**14, 1737000, "Luna")]);
"""
p'1 = v1
v'1 = Gm2(p2 - p1)/r**3_12
p'2 = v2
v'2 = Gm1(p1 - p2)/r**3_12

p1(0) = [0,0,0]
p2(0) = [362600., 0., 0.]
p'1(0) = [0., 0., 0.]
p'2(0) = [0., 1078.2, 0.]

dim = 5
W = [0, 0, 0, 0, 362600 0, 0, 0, 0, 0, 0, 1078.2, 0]
"""

def main():
    dt = 1./30 # 30 frames per second
    tol = 05e-10;
    sys.estimate = RungeKuttaFehlberg54(sys.ydot, 13, dt, tol);
                    #RungeKuttaFehlberg54(sys.ydot, sys.dim, dt, tol)];

    #F = lambda x

    #W, E = rkf54.safeStep(W);
    #rkf54.setStepLength(tEnd - W[0]);
    #W, E = rkf54.step(W);

    #print(W, E);
    #print("Total error: ", [abs(W[1] - y_1(1)), abs(W[2] - y_2(1))]);

    # The figure is set
    fig = plot.figure();
    axes = fig.add_subplot(111, aspect="equal", autoscale_on=False, xlim=(-sys.objects[1].position[0] - 50000, sys.objects[1].position[0] + 50000), ylim=(-sys.objects[1].position[0] - 50000, sys.objects[1].position[0] + 50000))
    line1, = axes.plot([], [], "o-b", lw=2); # Terra
    line2, = axes.plot([], [], "o-k", lw=2); # Luna
    time_text = axes.text(0.02, 0.95, "", transform=axes.transAxes);

    def init():
        #initialize animation
        line1.set_data([], []);
        line2.set_data([], []);
        time_text.set_text('');
        return line1, line2, time_text;

    def animate(i):
        #perform animation step
        global sys;
        sys.step(1);
        line1.set_data(*sys.objects[0].position[0:2]);
        line2.set_data(*sys.objects[1].position[0:2]);
        time_text.set_text('time = %.1f' % sys.time_elapsed());
        return line1, line2, time_text;

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