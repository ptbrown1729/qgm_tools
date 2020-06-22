function [H,V,xpts,ypts] = Hsho2D()
%Diagonalize 2D real space Hamiltonian
%H, the position space Hamiltonian is returned using a 'dictionary' basis.
%i.e. an eigenvector is made up of basis vectors in the following order...
%(x1,y1),(x1,y2),...(x1,yn),(x2,y1),...(xn,yn)
 
run('constants.m');
omega = 2*pi*3e3;
%define potential
V=@(x,y)0.5*mLi*omega^2*(x.^2+y.^2);
achar=sqrt(hbar/(2*mLi*omega));

%create x-grid
xstart = -5*achar;
xend = 5*achar;
dx = achar/10;
xpts = xstart:dx:xend;
Nx = length(xpts);
%create y-grid
ystart = -5.5*achar;
yend = 5.5*achar;
dy = achar/10;
ypts = ystart:dy:yend;
Ny = length(ypts);
%build position vector
pos = [kron(xpts',ones(length(ypts),1)),kron(ones(length(xpts),1),ypts')];
%(x1,y1),(x1,y2),...(x1,yn),(x2,y1),...
%Build Hamiltonian

Vmatr = sparse(1:Nx*Ny,1:Nx*Ny,V(pos(:,1),pos(:,2)));

DxForm = -2*sparse(1:Nx,1:Nx,ones(1,Nx))+sparse(1:(Nx-1),2:Nx,ones(1,Nx-1),Nx,Nx)+sparse(2:Nx,1:(Nx-1),ones(1,Nx-1),Nx,Nx);
DxSqrd = kron(DxForm,sparse(1:Ny,1:Ny,ones(1,Ny)))/dx^2;
DyForm = -2*sparse(1:Ny,1:Ny,ones(1,Ny))+sparse(1:(Ny-1),2:Ny,ones(1,Ny-1),Ny,Ny)+sparse(2:Ny,1:(Ny-1),ones(1,Ny-1),Ny,Ny);
DySqrd= kron(sparse(1:Nx,1:Nx,ones(1,Nx)),DyForm)/dy^2;

H = -0.5*hbar^2/(2*mLi)*(DxSqrd+DySqrd) + Vmatr; %()*

%diagonalize Hamiltonian
% [psivects,Es] = eigs(H,Neigs,'lm'); %'sm' is "smallest magnitude". 'lm' is largest magnitude.
% [Es,sortedIndices] = sort(diag(Es));
% psivects = psivects(:,sortedIndices);



end

