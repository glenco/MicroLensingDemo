# MicroLensingDemo

![On the left is the inverse magnification and on the right is an image of a lensed Gaussian source.  The convergence in stars is 1 here.](image.png)

This is a simple program to demonstrate how to do microlensing calculations with GLAMER.

There are more sophisticated things you can do with a `Grid` instead of a `GridMap` as used here 
in which case the grid can by dynamically refined the caustics or around the source 
images.  But if you have enough memory and enough cores the `GridMap` works just fine.

In your simulation you should make sure:

1) You have set the N_THREADS compiler flag to as many cores as you can 
spare (`cmake .. -DN_THREADS=20` in the build directory of GLAMER for example).  
This code is highly parallelized.

2) You have made a `GridMap` or `Grid` with a high enough initial resolution to 
resolve all the relavent images.

3) You don't move the source far enough that the images are no longer on 
the `Grid`.  You can always construct another grid centered on another 
point along the sources path.

Good luck!
