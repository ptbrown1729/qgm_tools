% TODO: this script seems broken ...

%Fundamental constants
c = 299792458; % m/s
h = 6.6261e-34; % Js
hbar = h/(2*pi); % Js
kb = 1.3806e-23; % J/K

%Lithium-6
mLi = 9.9883e-27; % kg
gamma = (2*pi) * 6e6; % Hz
lambda0 = 671e-9; 
f0 = c / lambda0;
w0 = (2*pi) * f0;

%Laser Params
lambda = 1070e-9;
f = c / lambda;
w = 2 * pi * f;
k = 2*pi /lambda;
Er = hbar^2 * k^2 / ( 2 * mLi );

%measured Bottom Beam vals
P = 2.25; %Watts % 2.25 W at 7V servo
wx0 = 104e-6; %41e-6;%108.7e-6;
zx = pi*wx0^2/lambda;
wy0 = 100;% 111.6e-6;
zy = pi*wy0^2/lambda;

%alpha = 0.00255; %Polarizability of 2s state for 1070nm in Hz/(W/m)^2

Int = 2*P/(pi*wx0*wy0);%*1e-10 %in MW/cm2
Vdip_Atom = 3*pi*c^2/(2*w0^3)*(gamma/(w0 - w) + gamma/(w0 + w))*Int;
%Vdip_Atom = h*alpha/2*Int; %From Full Calculation
Vdip_Mol = 2*Vdip_Atom;

%K = 3*pi*c^2/(2*w0^3)*(gamma/(w0 - w) + gamma/(w0 + w))*(2*P/pi)*(4/mLi);
%Trap depth and trapping frequencies
%I(x,y,z) = (2P/pi wx0 wy0)*(1 - 2 x^2/wx0^2 - 2 y^2/wy0^2 - 0.5 z^2*(1/zx^2+1/zy^2)
%assuming our beam is not astigmatic.
%Depth depends on atoms vs. molecules, but trapping freq. doesn't.
DepthMol = Vdip_Mol/kb;
DepthAtoms = Vdip_Atom/kb;
fx0 = sqrt(4*Vdip_Atom/(mLi*wx0^2))/(2*pi);
fy0 = sqrt(4*Vdip_Atom/(mLi*wy0^2))/(2*pi);
fz0 = sqrt(2*Vdip_Atom*0.5*(1/zx^2+1/zy^2)/mLi)/(2*pi);

PCakeDeConfX = 100.4; %Guess for pancake radial deconfinement (Hz) 
PCakeDeConfY = 77.33; %at maximum axial trap freq (Amp=3.2V, Freq=7V)

fx0eff = sqrt(fx0^2-PCakeDeConfX^2);
fy0eff = sqrt(fy0^2-PCakeDeConfY^2);

disp(cat(2,num2cell([DepthMol,DepthAtoms,fx0,fy0,fz0,fx0eff,fy0eff])',{'Depth Molecules','Depth Atoms','fx0','fy0','fz0','fx0eff','fy0eff'}'))

%if beam is astigmatic, modify trapping frequencies.
fx=@(z)omegax0./(1+(z/zy).^2).^(1/4);
fy=@(z)omegay0./(1+(z/zx).^2).^(1/4);