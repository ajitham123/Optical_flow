clc; clear all; close all;

% load images:
Im1 = imread('images/frames99.jpg');
Im2 = imread('images/frames100.jpg');

% determine the individual flow vectors:
graphics = true;
im_scale = true;
[point, flow_mag, angle] = get_optical_flow_edges(Im1, Im2, graphics, im_scale);

% find FOE
FOE = find_FOE(point, flow_mag, angle,graphics)

