HOW TO USE:

To find and plot velocity inhibition during trajectory intersection events:
1. Run "vidTrack.m" and select tif image stack.
2. Run "findCrossings.m" to find events where trajectories intersect.
3. Run "intersectionPlot.m" to calculate the velocity reduction at intersection events and plot the reduction in velocity vs. dt.

To generate SI video:
1. Run "vidTrackSIVid.m" and select "stacks/4.1 6x 15fps_5x speed.tif" to find trajectories seen in the SI video.
2. Run "findCrossingsSIVid.m" to identify intersection events. 
3. Run "videoGeneratorSIVid.m" to generate frames for video with particle trajectories overlayed on top of the raw video.