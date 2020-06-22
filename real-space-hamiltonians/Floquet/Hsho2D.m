function [H] = Hsho2D(t)
%Hamiltonianian for 2D SHO with rotating anisotropy
%in number basis |nx ny>
run('constants.m')
%states from nx=0 to nx=Nx
Nx=10;
Ny=11;
omega=2*pi*3e3;
epsilon=2*pi*3e2;
modW=2*pi*6e3;

H0=hbar*omega*(kron(diag(0:1:Nx),eye(Ny+1))+kron(eye(Nx+1),diag(0:1:Ny)));
HyconnDiag=hbar*epsilon/2*cos(modW*t)*kron(diag(ones(Nx+1,1)),diag(0:1:Ny));
HyconnOffDiag=hbar*epsilon/2*cos(modW*t)*kron(diag(ones(Nx+1,1)),diag(sqrt(1:Ny-1).*sqrt(2:Ny),2));

HxconnDiag=hbar*epsilon/2*sin(modW*t)*kron(diag(0:1:Nx),diag(ones(Ny+1,1)));
HxconnOffDiag=hbar*epsilon/2*sin(modW*t)*kron(diag(sqrt(1:Nx-1).*sqrt(2:Nx),2),diag(ones(Ny+1,1)));

H=H0+HyconnDiag+HyconnOffDiag+HyconnOffDiag'+HxconnDiag+HxconnOffDiag+HxconnOffDiag';

end

