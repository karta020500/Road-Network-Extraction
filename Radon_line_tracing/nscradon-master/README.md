# Detect lines on needle shaped crystals images

This repository contains several MATLAB functions and a demo script we used to detect needle shaped 
crystals on images for monitoring their growth. The results and description of the method can be 
found in [1]. Please reference this publication if you use this code for your research.

The method is based on MATLAB's 
[Radon transform based line detector](https://se.mathworks.com/help/images/detect-lines-using-the-radon-transform.html)
with several improvements in order to increase its perfromance and quality for the images with large 
amount of line segments. The main imrpovements are following:

1. Splitting image into segments (using `bwlabel()`) and processing the segments separately
2. Using Gaussian blur to remove a "baseline" for better peak detection in Rho/Tau space
3. Possibility for fine tuning of the algorithm 

The main function is `getlinesforbw.m`, it takes a binary image as a main argument
and returns a structure array with detected lines. You can also provide a structure 
with parameters as a second argument, see `demo.m` for example and more details.

Notice, that:

1. Changing some of the parameters (e.g. `theta_step`) increases computational time.
2. Optimization of the paremeters using e.g. DoE approach is a good idea.

The scripts require Image Processing Toolbox.

## References
1. [Image Analytical Approach for Needle-Shaped Crystal Counting and Length Estimation](https://pubs.acs.org/doi/10.1021/acs.cgd.5b00720). 
Wu, Jian X.; Kucheryavskiy, Sergey V.; Jensen, Linda G.; Rades, Thomas; Müllertz, Anette; 
Rantanen, Jukka. Crystal Growth & Design, Bind 15, Nr. 10, 2015, s. 4876-4885.
