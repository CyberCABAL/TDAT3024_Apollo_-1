from fy import CelestialObject
from orbitalsys import System
from RungeKuttaFehlberg import RungeKuttaFehlberg54

import numpy as np
import time
import matplotlib.pyplot as plot
import matplotlib.animation as animation
import matplotlib.patches as patch

Terra = CelestialObject([0, 0, 0], [0, 0, 0], 59722 * 10**20, 6378100, "Terra");
Luna = CelestialObject([362600.0, 0, 0], [0, 1078.2, 0], 734767309 * 10**14, 1737000, "Luna");
sys = System([Terra, Luna]);

def main():
    dt = 1./30 # 30 frames per second

    W = np.array([0, 1, 0]);
    #h = 0.25;
    tol = 05e-14;
    #tEnd = 100.0;
    #F = lambda x
    #rkf54 = RungeKuttaFehlberg54(F, 3, h, tol);

    #while(W[0] < tEnd):
    #    W, E = rkf54.safeStep(W);
        
    #rkf54.setStepLength(tEnd - W[0]);
    #W, E = rkf54.step(W);

    #print(W, E);
    
    #print("Total error: ", [abs(W[1] - y_1(1)), abs(W[2] - y_2(1))]);

    # The figure is set
    fig = plot.figure();
    axes = fig.add_subplot(111, aspect="equal", autoscale_on=False, xlim=(-Luna.position[0] - 50000, Luna.position[0] + 50000), ylim=(-Luna.position[0] - 50000, Luna.position[0] + 50000))
    
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
        global dt, rkf54, Terra, Luna, sys;
        #for j in range(10):
            #rkf54.safeStep(dt);
        line1.set_data(*Terra.position[0:2]);
        line2.set_data(*Luna.position[0:2]);
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
