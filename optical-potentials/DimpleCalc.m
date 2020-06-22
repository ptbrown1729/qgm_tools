%Lightsheet depth and trapping frequencies. Waists given here are based off
%of measurged values.

%Fundamental constants
c = 3*10^8;
h = 6.63*10^(-34);
hbar = h/(2*pi);
kb = 1.38*10^(-23);

%Lithium-6
mLi = 6*1.67*10^(-27);
gamma = 2*pi*6*10^6;
lambda0 = 671*10^(-9);
f0 = c/lambda0;
w0 = 2*pi*f0;

%LS Laser
lambda = 832*10^(-9);
f = c/lambda;
w = 2*pi*f;
k = 2*pi/lambda;
Er = hbar^2*k^2/(2*mLi);

%measured LS vals
P = 28e-3;
wx0 = 12.5e-6;
wy0 = wx0;%49e-6;%125*10^(-6);
zx = pi*wx0^2/lambda;
zy = pi*wy0^2/lambda;

Int = 2*P/(pi*wx0*wy0);%*1e-10 %in MW/cm2
Vdip_Atom = 3*pi*c^2/(2*w0^3)*(gamma/(w0 - w) + gamma/(w0 + w))*Int;
Vdip_Mol = 2*Vdip_Atom;

%K = 3*pi*c^2/(2*w0^3)*(gamma/(w0 - w) + gamma/(w0 + w))*(2*P/pi)*(4/mLi);
%Trap depth and trapping frequencies
%I(x,y,z) = (2P/pi wx0 wy0)*(1 - 2 x^2/wx0^2 - 2 y^2/wy0^2 - 0.5 z^2*(1/zx^2+1/zy^2)
%assuming our beam is not astigmatic.
%Depth depends on atoms vs. molecules, but trapping freq. doesn't.
DepthAtoms = Vdip_Atom/kb;
fx0 = sqrt(4*Vdip_Atom/(mLi*wx0^2))/(2*pi);
fy0 = sqrt(4*Vdip_Atom/(mLi*wy0^2))/(2*pi);
fz0 = sqrt(2*Vdip_Atom*0.5*(1/zx^2+1/zy^2)/mLi)/(2*pi);

fprintf('Depth Atoms = %0.1f KHz \n',Vdip_Atom/hbar/10^3);
fprintf('Depth Atoms = %0.1f uK \n',Vdip_Atom/kb/1e-6);
fprintf('Omega_x = (2pi) %0.2f KHz \n',fx0/1e3);
fprintf('Omega_y = (2pi) %0.2f KHz \n',fy0/1e3);
fprintf('Omega_z = (2pi) %0.2f KHz \n\n',fz0/1e3);

%disp(cat(2,num2cell([DepthMol,DepthAtoms,fx0,fy0,fz0])',{'Depth Molecules','Depth Atoms','fx0','fy0','fz0'}'))

%if beam is astigmatic, modify trapping frequencies.
fx=@(z)omegax0./(1+(z/zy).^2).^(1/4);
fy=@(z)omegay0./(1+(z/zx).^2).^(1/4);