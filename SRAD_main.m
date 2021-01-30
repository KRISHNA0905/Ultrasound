clc;
clear all;
close all;

sigma = 0.4;
I = phantom(256);

I = double(20*(I+0.5));
N = (randn(size(I)).*sigma);
img1 = I+(I).*N;
%img1 = imresize(img1,0.5);
imagesc(img1);
colormap gray;
[out_img rec]=SRAD(img1,1000,0.2);
figure;
imagesc(out_img);
colormap gray;
loss=SNR(img1,out_img)