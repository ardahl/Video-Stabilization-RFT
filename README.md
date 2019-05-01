# Video Stabilization #
This was a class project I did to try to do video stabilization using the [Robust Feature Trajectories](https://www.csie.ntu.edu.tw/~cyy/publications/papers/Lee2009VSR.pdf) method. It also includes a method for simple video stabilization by matching transformations frame-by-frame and smoothing them. RTF is in the StabilizeVideoRTF.m and the simple version is in StabilizeVideoSimple.m.

A tangental part of this project was to do an analysis on the error of smoothing transformation matrices. The smoothing was done by decomposing the transform into a translation, rotation, and scale, averaging those individually, and reforming the transformation. This was compared to averaging the components of the transformations directly to see where the error manifested.
