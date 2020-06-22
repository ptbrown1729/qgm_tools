% Calculate trapping frequencies, depth, etc. for a single beam

% Fundamental constants
c = 299792458; % m/s
h = 6.6261e-34; % Js
hbar = h/(2*pi); % Js
kb = 1.3806e-23; % J/K

% Lithium-6
mLi = 9.9883e-27; % kg
gamma = (2*pi) * 6e6; % Hz
lambda0 = 671e-9; 
f0 = c / lambda0;
w0 = (2*pi) * f0;

% Laser Params
lambda = 1070e-9;
f = c / lambda;
w = 2 * pi * f;
k = 2*pi /lambda;
Er = hbar^2 * k^2 / ( 2 * mLi );

%Beam Params
P = 2.25; % W
% waists (2 * sigma)
wx0 = 105e-6; %68e-6;
wy0 = 99e-6; %wx0;
% rayleigh ranges
zx = pi * wx0^2 / lambda;
zy = pi * wy0^2 / lambda;

% Potential
peak_intensity = 2 * P / (pi*wx0*wy0); %*1e-10 %in MW/cm^2
Vdip_Atom = 3*pi*c^2 / (2*w0^3) * (gamma/(w0 - w) + gamma/(w0 + w)) * peak_intensity;
Vdip_Mol = 2 * Vdip_Atom;

%Trap depth and trapping frequencies
%I(x,y,z) = (2P/pi wx0 wy0)*(1 - 2 x^2/wx0^2 - 2 y^2/wy0^2 - 0.5 z^2*(1/zx^2+1/zy^2)
%assuming our beam is not astigmatic.

% trap depth depends on atoms vs. molecules, but trapping frequency doesn't.
fx0 = sqrt( 4 * abs(Vdip_Atom) / (mLi*wx0^2) ) / (2*pi);
fy0 = sqrt( 4 * abs(Vdip_Atom) / (mLi*wy0^2) ) / (2*pi);
fz0 = sqrt( 2 * abs(Vdip_Atom) * 0.5 * (1/zx^2+1/zy^2) / mLi ) / (2*pi);

% scattering rate
% Gamma = (2pi) frq
scatt_rate = (1/hbar) * 3 * pi * c^2 / (2*w0^3) * (w/w0)^3 * (gamma/(w0 - w) + gamma/(w0 + w))^2 * peak_intensity;


fprintf('Depth Atoms = %0.1f uK  = %0.0f KHz\n', Vdip_Atom / kb / 1e-6, Vdip_Atom / h / 10^3);
fprintf('Omega_x = (2pi) %0.1f Hz \n', fx0);
fprintf('Omega_y = (2pi) %0.1f Hz \n', fy0);
fprintf('Omega_z = (2pi) %0.1f Hz \n', fz0);
fprintf('Scattering Rate = (2pi) %0.0f mHz \n', scatt_rate/(2*pi)*1e3);
fprintf('Scattering time = %0.0f s \n\n', 1/(scatt_rate/(2*pi)));
