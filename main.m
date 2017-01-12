% ---------------------------------------------------------------------------
%% Initialize
clear all
clear
clc
close all
addpath('./InterInner');

%% Demosaicking
i=19;

filename=['kodim' num2str(i) '.png']
img = imread(filename);

img_InnerInter = demosaicking_InnerInter(img);
imshow(uint8(img_InnerInter));