function [H,V,xpts] = HdoubleWell(P)
%[H,V,xpts] = HdoubleWell(P)
run('constants.m')

%define potential

% P = 90;
% w1 = 120e-6;
omega0 = c*2*pi/671e-9;
omega = c*2*pi/1064e-9;
gamma = 2*pi*6e6;
alpha = ((3*pi*c^2)/(2*omega0^3))*(gamma/(omega0 - omega) + gamma/(omega0 + omega));
w1 = 120e-6;
w2 = 120e-6;
Pow1 = 2*1e-3;
Pow2 = 1e-3;
Cx1 = -0.5*w1;
Cx2 = 0.5*(w1+w2);
V=@(x) -alpha*(2*Pow1/(pi*w1^2))*exp(-2*(x-Cx1).^2/w1^2)-alpha*(2*Pow2/(pi*w2^2))*exp(-2*(x-Cx2).^2/w2^2)+alpha*(2*Pow1/(pi*w1^2))+alpha*(2*Pow2/(pi*w2^2));

Depth1 = alpha*(2*Pow1/(pi*w1^2));
OmegaSHO1 = sqrt(4*Depth1/(mLi*w1^2));
achar1 = sqrt(hbar/(2*mLi*OmegaSHO1));
Depth2 = alpha*(2*Pow2/(pi*w1^2));
OmegaSHO2 = sqrt(4*Depth2/(mLi*w2^2));
achar2 = sqrt(hbar/(2*mLi*OmegaSHO2));

achar = min(achar1,achar2);
wmax = max(w1,w2);
cmin = min(Cx1,Cx2);
cmax = max(Cx1,Cx2);

%initialize grid
xstart = cmin-2*wmax;
xend = cmax+2*wmax;
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

