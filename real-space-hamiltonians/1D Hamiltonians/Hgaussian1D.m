function [H,V,xpts,OmegaSHO] = Hgaussian1D(P,w0)
run('constants.m')

%define potential

% P = 90;
% w0 = 120e-6;
omega0 = c*2*pi/671e-9;
omega = c*2*pi/1064e-9;
gamma = 2*pi*6e6;
alpha = ((3*pi*c^2)/(2*omega0^3))*(gamma/(omega0 - omega) + gamma/(omega0 + omega));
V=@(x) -1*alpha*(2*P/(pi*w0^2))*exp(-2*x.^2/w0^2)+alpha*(2*P/(pi*w0^2));

Depth = alpha*(2*P/(pi*w0^2));
OmegaSHO = sqrt(4*Depth/(mLi*w0^2));
achar = sqrt(hbar/(2*mLi*OmegaSHO));

%initialize grid
xstart = -3*w0;
xend = 3*w0;
xspace = achar/10;
xpts = xstart:xspace:xend;
Nx = length(xpts);

%define Hamiltonian
Potential = sparse(1:Nx,1:Nx,V(xpts));
%Potential = diag(V(xpts));
%dx = diag(-1/(2*xspace)*ones(1,Nx-1),-1) + diag(1/(2*xspace)*ones(1,Nx-1),1);
dxSqr =(1/xspace^2)*(-2*sparse(1:Nx,1:Nx,ones(1,Nx))+sparse(1:Nx-1,2:Nx,ones(1,Nx-1),Nx,Nx)+sparse(2:Nx,1:Nx-1,ones(1,Nx-1),Nx,Nx));
H = -hbar^2/(2*mLi)*dxSqr + Potential;


end