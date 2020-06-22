function [H,V,xpts] = HlatticeWithTrap(omega,LattDepth,LattSpacing)
%function [H] = Hsho1D(omega,LattDepth,LattSpacing)
%V=@(x) 0.5*mLi*omega^2*x.^2 + LattDepth/2*(1 - cos(k*x))

run('constants.m')

%define potential
k = 2*pi/LattSpacing;
V=@(x) 0.5*mLi*omega^2*x.^2 + LattDepth/2*(1 - cos(k*x));
achar = sqrt(hbar/(2*mLi*omega));

%initialize grid
xstart = -100*achar;
xend = 100*achar;
%xspace = achar/10;
xspace = LattSpacing/20;
%Nx = round((xend-xstart)/xspace+1);
xpts = xstart:xspace:xend;
Nx = length(xpts);

%define Hamiltonian
Potential = sparse(1:Nx,1:Nx,V(xpts));
%Potential = diag(V(xpts));
%dx = diag(-1/(2*xspace)*ones(1,Nx-1),-1) + diag(1/(2*xspace)*ones(1,Nx-1),1);
dxSqr =(1/xspace^2)*(-2*sparse(1:Nx,1:Nx,ones(1,Nx))+sparse(1:Nx-1,2:Nx,ones(1,Nx-1),Nx,Nx)+sparse(2:Nx,1:Nx-1,ones(1,Nx-1),Nx,Nx));
H = -hbar^2/(2*mLi)*dxSqr + Potential;


end


