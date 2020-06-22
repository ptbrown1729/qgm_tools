%06/07/2016
%Model for non-interacting Fermi gas in a lattice. Using tight-binding
%limit to solve for 1D density of states, and then get atomic density in
%LDA using semi-classical approach. Need U<4t or so for this approximation
%to become reasonable.

h = 6.626e-34;
hbar = h/(2*pi);
mLi = 9.9883e-27;
omega = 2*pi*140;

%for this to work appropriately, must have Ef<4t. Otherwise not single band
%and approximation breaks.
a = 532e-9*sqrt(2);
t = h*1e3;
Ef = h*30e3;
OmegaZ = 2*pi*20e3;

V = @(r) 0.5*mLi*omega^2*r.^2;
if Ef>4*t
    Ef_eff1 = 4*t;
else
    Ef_eff1=Ef;
end
if Ef>(4*t+hbar*OmegaZ)
    Ef_eff2 = 4*t+hbar*OmegaZ;
else
    Ef_eff2 = Ef;
end

n1stAxialState = @(r) acos(1-(Ef_eff1-V(r))/(2*t))*4/pi/a;
n2ndAxialState = @(r) acos(1-(Ef_eff2-hbar*OmegaZ-V(r))/(2*t))*4/pi/a*heaviside(Ef-hbar*OmegaZ);
Rpts = linspace(0,120e-6,200);

N1st = n1stAxialState(Rpts);
N1st(imag(N1st)>1e-30)=0;
N2nd = n2ndAxialState(Rpts);
N2nd(imag(N2nd)>1e-30)=0; 
Nprofile = N1st+N2nd; 

plot(Rpts,Nprofile)
hold on;