close all
clc

 
video = VideoWriter('VideoHINOCRE.mp4','MPEG-4'); %create the video object
video.FrameRate = 20; % How many frames per second.
open(video); %open the file for writing
 
for ii=1:2728 %where N is the number of images 
    
  I = imread(['Video2\' num2str(ii,'%04.f') '.png']); %read the next image 
 
  writeVideo(video,I); %write the image to file
%     delete(I) 
end
clearvars -except video
close(video); %close the file



  
function [imageLeft,imageRight] = stereoRead(path,CamRead,k,stereoParams)
% read and rectify with stereoParams
imageLeft  = imread([path '\cam0\data\',num2str(CamRead(k),'%d'),'.png']);
imageRight = imread([path '\cam1\data\',num2str(CamRead(k),'%d'),'.png']);
[imageLeft,imageRight] = rectifyStereoImages(imageLeft,imageRight,stereoParams);
end
