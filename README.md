# Optical_flow
MAVLab, TU Delft

Ajith Anil Meera, Guido de Croon

Abstract—The contemporary techniques employed to estimate the optical flow parameters between a pair of consecutive images utilize the entire image pixels for image correspondences using non linear optimization techniques, which makes the algorithm computationally heavy, especially for real time applications with on-board computing, like on drones for example. This paper presents an edge based algorithm to evaluate the flow parameters in a way that is sparse and computationally cheap - without having to match every pixel in the image - by evaluating those edge pixels in the image with highest confidence for a flow to have occurred along the direction of the tangent to the edge at that point. The algorithm handles the noise in the estimated flow vectors from disrupting and destabilizing the solution by iteratively removing them. The algorithm was validated by artificially generating an image database with known flow parameters and comparing the results of the algorithm to it. The results of the algorithm were found to be better than the Lukas Kanade method applied on Harris corners at the expense of higher computational time.

MATLAB codes:

demo.m - Run this

get_optical_flow_point.m - Function that finds the optical flow at the edges in the direction of gradient

find_foe.m - Function to find the FOE

remove_outlier.m - Function to remove outlier tangent lines with zero flow magnitude

Lucas_Kanade.m - Code to generate histogram and heat map for FOE error using LK method on Harris corners

C codes:

main.cpp inside the c_optical_flow folder

Documentation:

report.pdf - Internship report

