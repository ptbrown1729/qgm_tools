% created 2018/12/04 by pb
% Consider a lattice formed by two beams intersecting at angle 2 * theta


%Fundamental constants
c = 299792458; % m/s
h = 6.6261e-34; % Js
hbar = h/(2*pi); % Js
kb = 1.3806e-23; % J/K

%Lithium-6
mLi = 9.9883e-27; % kg
gamma = 2*pi * 6e6; % Hz
lambda0 = 671e-9; % m
f0 = c / lambda0; % Hz
w0 = 2 * pi * f0; % Hz

%Laser Params
lambda = 532e-9; % m
f = c/lambda; 
w = 2 * pi * f;
k = 2 * pi / lambda;

% angle
half_angle = 10 * pi/180; % radians
waist_x = 190-6; % m, these from measurements on 12/5/2018
waist_y = 30e-6; % m
zx = pi * waist_x ^ 2 / lambda;
zy = pi * waist_y ^ 2 / lambda;

% lattice spacing
latt_spacing = lambda / ( 2 * sin(half_angle) );

