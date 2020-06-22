function [H] = H2D(t)
%H2D Summary of this function goes here
%   Detailed explanation goes here

%problem with this seems to be that exponentiating sparse matrices is not a
%good operation.

run('constants.m')

trapFrq=2*pi*3e3;
epsilon=2*pi*3e2;
omega=0;
V=@(x,y)hbar*0.5*((trapFrq+epsilon*cos(omega*t))^2*x.^2+(trapFrq+epsilon*sin(omega*t))^2*y.^2);
achar=sqrt(hbar/(2*mLi*trapFrq));

%create x-grid
xstart = -5*achar;
xend = 5*achar;
dx = achar/10;
xpts = xstart:dx:xend;
Nx = length(xpts);
%create y-grid
ystart = -5*achar;
yend = 5*achar;
dy = achar/10;
ypts = ystart:dy:yend;
Ny = length(ypts);
%build position vector
pos = [kron(xpts',ones(length(ypts),1)),kron(ones(length(xpts),1),ypts')];
%(x1,y1),(x1,y2),...(x1,yn),(x2,y1),...
%Build Hamiltonian

Vmatr = sparse(1:Nx*Ny,1:Nx*Ny,V(pos(:,1),pos(:,2)));

DxForm = sparse(1:Nx,1:Nx,ones(1,Nx))*2+sparse(1:Nx-1,2:Nx,ones(1,Nx-1),Nx,Nx)+sparse(2:Nx,1:Nx-1,ones(1,Nx-1),Nx,Nx);
DxSqrd = kron(DxForm,sparse(1:Ny,1:Ny,ones(1,Ny)))/dx^2;
DyForm = sparse(1:Ny,1:Ny,ones(1,Ny))*2+sparse(1:Ny-1,2:Ny,ones(1,Ny-1),Ny,Ny)+sparse(2:Ny,1:Ny-1,ones(1,Ny-1),Ny,Ny);
DySqrd= kron(sparse(1:Nx,1:Nx,ones(1,Nx)),DyForm)/dy^2;

H = -0.5*(hbar^2/(2*mLi))*(DxSqrd+DySqrd) + Vmatr;


end

