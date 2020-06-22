function [H,V,xpts] = Hcoulomb1D()
run('constants.m')

%define potential

%need to exclude x=0 case...
V=@(x) (-e./x.^2);


%initialize grid
xstart = -10.000001*abohr;
xend = 10.000001*abohr;
xspace = abohr/10;
Nx = round((xend-xstart)/xspace+1);
xpts = xstart:xspace:xend;

%define Hamiltonian
Potential = sparse(1:Nx,1:Nx,V(xpts));
%Potential = diag(V(xpts));
%dx = diag(-1/(2*xspace)*ones(1,Nx-1),-1) + diag(1/(2*xspace)*ones(1,Nx-1),1);
dxSqr =(1/xspace^2)*(-2*sparse(1:Nx,1:Nx,ones(1,Nx))+sparse(1:Nx-1,2:Nx,ones(1,Nx-1),Nx,Nx)+sparse(2:Nx,1:Nx-1,ones(1,Nx-1),Nx,Nx));
H = -hbar^2/(2*mLi)*dxSqr + Potential;


end


