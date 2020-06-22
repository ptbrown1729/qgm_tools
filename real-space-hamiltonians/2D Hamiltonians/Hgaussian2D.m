function [H,V,xpts,ypts] = Hgaussian2D(P)
%[H,V,xpts,ypts] = Hgaussian2D(P)
%P = [Cx,Cy,Wx,Wy,Pow]
Cx = P(1); Cy = P(2);
Wx = P(3); Wy = P(4);
Pow = P(5);

run('..\constants.m')

omega0 = c*2*pi/671e-9;
omega = c*2*pi/1064e-9;
gamma = 2*pi*6e6;
alpha = ((3*pi*c^2)/(2*omega0^3))*(gamma/(omega0 - omega) + gamma/(omega0 + omega));
V=@(x,y) -alpha*(2*Pow/(pi*Wx*Wy))*exp(-2*(x-Cx).^2/Wx^2-2*(y-Cy).^2/Wy^2)+alpha*(2*Pow/(pi*Wx*Wy));

Depth = alpha*(2*Pow/(pi*Wx*Wy));
OmegaSHOX = sqrt(4*Depth/(mLi*Wx^2));
OmegaSHOY = sqrt(4*Depth/(mLi*Wy^2));
acharX = sqrt(hbar/(2*mLi*OmegaSHOX));
acharY = sqrt(hbar/(2*mLi*OmegaSHOY));

%create x-grid
xstart = -Wx;
xend = Wx;
dx = acharX/10;
xpts = xstart:dx:xend;
Nx = length(xpts);
%create y-grid
ystart = -Wy;
yend = Wy;
dy = acharY/10;
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

H = -0.5*hbar^2/(2*mLi)*(DxSqrd+DySqrd) + Vmatr;



end