close all;
clc;

y1d = -5 : 0.1 : 5;
x1d = -30 : 0.5 : 30;
[yy, xx] = meshgrid(x1d, y1d);

% figure
% % imagesc(x)

phi1 = 4*pi/180;
% phi2 = 4*pi/180;
phi2 = pi*0.0;

lsdepth = 3; % uK
greenDepth = 0.4; % uK

% lightsheet parameters
sigmax = 5.1;%um
sigmay = 49;%um
% ls = sin(x*cos(phi1)+y*sin(phi1)).^2;
% This looks wrong bc using should be using waists but using sigmas instead
ls = lsdepth * exp( -xx .^2 / 2 / sigmax^2 - yy.^2 / 2 / sigmay^2);

greenspacing = 2; %um
green = greenDepth * sin( 1/greenspacing * pi * ( xx * cos(phi1) + yy * sin(phi1)) + phi2 ) .^ 2;

combined = ls - green;

figure;

imagesc(x1d,y1d,ls);
axis equal;

figure;
imagesc(x1d,y1d,green);
axis equal;

figure;
imagesc(x1d,y1d,combined);
axis equal;

ref = max(max(combined))
atomfill = combined>ref*0.93;

figure
subplot(2,2,1);
imagesc(x1d,y1d,atomfill);
axis equal;

subplot(2,2,2);
plot(x1d,sum(atomfill,1))
