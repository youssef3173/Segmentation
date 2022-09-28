clear; close all; clc; 

img = imread('materiau.bmp');
[row, col] = size(img); 
figure(1), imshow(img); 
figure(2),  imhist(img);
imgTreshold = grayslice(img,14); 
figure(3), imshow(imgTreshold, colormap(jet(3)))

beta = 1; 
Mu = [10 130 190]; 
Sigma = [15 15 15]; 
EtiqInit = zeros(row, col); 
NbClasses = 3; 
EtiqResult=icm(img, EtiqInit, NbClasses, Mu, Sigma, beta);
%EtiqResult2=recuitSim(img, EtiqInit, NbClasses, Mu, Sigma, beta); 
%figure(4), imagesc(EtiqResult); 
figure(4), imshow(EtiqResult2, jet(4))