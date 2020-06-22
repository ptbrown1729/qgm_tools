function [H,xpts] = Hlattice(V,a) %Hlattice(V,a,Nbrill);
%Returns real-space lattice Hamiltonian.

run('constants.m')
%define potential
% lambda = 1064e-9;
% k=2*pi/lambda;
% Er=hbar^2*k^2/(2*mLi);
% V=@(x) 0.2*Er*cos(2*k*x);
%repeats with period lambda/2

Nbrill = 201;
%initialize grid
%dx should be a fraction of a.
dx = a/50; 
% xstart = 0;
% xend = Nbrill*a-dx;
xstart = -Nbrill*a/2;
xend = Nbrill*a/2-dx;
%Nx = round((xend-xstart)/dx+1);
xpts = xstart:dx:xend;
Nx = length(xpts);

%define Hamiltonian
Potential = sparse(1:Nx,1:Nx,V(xpts));
%Potential = diag(V(xpts));
%Periodic second derivative operator
dxSqrPeriodic =(1/dx^2)*(-2*sparse(1:Nx,1:Nx,ones(1,Nx))+sparse(1:Nx-1,2:Nx,ones(1,Nx-1),Nx,Nx)+sparse(2:Nx,1:Nx-1,ones(1,Nx-1),Nx,Nx)+sparse(1,Nx,1,Nx,Nx)+sparse(Nx,1,1,Nx,Nx));
H = -hbar^2/(2*mLi)*dxSqrPeriodic + Potential;

end

