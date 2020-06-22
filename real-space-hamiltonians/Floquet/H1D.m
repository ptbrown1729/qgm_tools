function [H] = H1D(T)

run('constants.m')
lambda=1064e-9;
omega=2*pi*c/lambda;
kr=2*pi/lambda;
b=lambda/10; %maximum lattice displacement
omegaFl=(2*pi)*3e3;

I=2*0.1/(pi*(100e-6)^2);
V0=(3*pi*c^2)/(2*omega_D2^3)*(lambda_D2/(omega_D2-omega)+lambda_D2/(omega_D2+omega))*I;
V=@(x) V0*cos(kr*x).^2;

%initialize grid
xstart = -3*lambda;
xend = 3*lambda;
dx = lambda/100;
N = (xend-xstart)/dx+1;
xpts = xstart:dx:xend;

%define Hamiltonian


Potential = diag(V(xpts));
DX = 1/dx*diag(-0.5*ones(1,N-1),-1) + diag(0.5*ones(1,N-1),1);
DXSqr = 1/dx^2*(diag(ones(1,N-1),-1) + diag(ones(1,N-1),1) - 2*diag(ones(1,N)));
H = -0.5*hbar^2/(2*mLi)*DXSqr + Potential -b*omegaFl*sin(omegaFl*T)*(-1i*DX);


end

