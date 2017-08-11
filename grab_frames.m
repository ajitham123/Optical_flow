Folder = 'C:\Users\Ajith A M\OneDrive\Windows desktop backup\MAVlab\frames';
clc;
clear all;
close all;
tic;
vid=VideoReader('VID_20170801_162240862.mp4');
numFrames = vid.NumberOfFrames;
n=numFrames;
for i = 1:1:n
  frames = read(vid,i);
  imwrite(frames,['C:\Users\Ajith A M\OneDrive\Windows desktop backup\MAVlab\Optical_flow\frames\frames', int2str(i), '.jpg']);
end 
