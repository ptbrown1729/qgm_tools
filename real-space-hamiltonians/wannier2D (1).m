%constants
c=3e8;
hbar = 1.054e-34;
m = 9.9883e-27;
kb = 1.38e-23;

%create potential
lambda = 1064e-9;
kl = 2*pi/lambda;
omega = c*kl;
lambda0 = 671e-9;
omega0 = c*2*pi/lambda0;
gamma = 2*pi*6e6;

P=0.1;
w0=100e-6;
zr = pi*w0^2/lambda;
w = @(x) w0*sqrt(1+(x/zr).^2);
V=@(x,y)(3*pi*c^2/(2*omega0^3)*(gamma/(omega0-omega)-gamma/(omega0+omega)))*(2*P./(pi*w(x).^2)).*exp(-2*y.^2./w(x).^2);

%create x-grid
xstart = -2*zr;
xend = 2*zr;
dx = zr/100;
xpts = xstart:dx:xend;
Nx = length(xpts);
%create y-grid
ystart = -2*w0;
yend = 2*w0;
dy = w0/50;
ypts = ystart:dy:yend;
Ny = length(ypts);
%build position vector
pos = [kron(xpts',ones(length(ypts),1)),kron(ones(length(xpts),1),ypts')];
%(x1,y1),(x1,y2),...(x1,yn),(x2,y1),...
%Build Hamiltonian
well = @(x,y) 0.5*(x.^2+y.^2-1).*(abs(x.^2+y.^2)<1);

Vmatr = sparse(1:Nx*Ny,1:Nx*Ny,V(pos(:,1),pos(:,2)));

DxForm = sparse(1:Nx,1:Nx,ones(1,Nx))*2+sparse(1:Nx-1,2:Nx,ones(1,Nx-1),Nx,Nx)+sparse(2:Nx,1:Nx-1,ones(1,Nx-1),Nx,Nx);
DxSqrd = kron(DxForm,sparse(1:Ny,1:Ny,ones(1,Ny)))/dx^2;
DyForm = sparse(1:Ny,1:Ny,ones(1,Ny))*2+sparse(1:Ny-1,2:Ny,ones(1,Ny-1),Ny,Ny)+sparse(2:Ny,1:Ny-1,ones(1,Ny-1),Ny,Ny);
DySqrd= kron(sparse(1:Nx,1:Nx,ones(1,Nx)),DyForm)/dy^2;

%DxForm = sparse(diag(2*ones(Ny,1))+diag(ones(Ny-1,1),1)+diag(ones(Ny-1,1),-1));
%DxSqrd = sparse(kron(DxForm,diag(ones(Nx,1))))/dx^2;
%DyForm = sparse(diag(2*ones(Nx,1))+diag(ones(Nx-1,1),1)+diag(ones(Nx-1,1),-1));
%DySqrd = sparse(kron(diag(ones(Ny,1)),DyForm))/dy^2;

H = -0.5*(hbar^2/(2*m))*(DxSqrd+DySqrd) + Vmatr;
%diagonalize Hamiltonian
Neigvects = 5;
[psivects,E] = eigs(H,Neigvects,'sm'); %note 'sm' is "smallest magnitude"
[E,sortedIndices] = sort(diag(E));
psivects = psivects(:,sortedIndices);

subplot(3,3,1)
imagesc(reshape(V(pos(:,1),pos(:,2)),Ny,Nx));
colorbar;
axis equal;

subplot(3,3,2)
plot([0,1],[E,E],'b')

for ii=1:5
subplot(3,3,ii+2)
imagesc(reshape(psivects(:,ii),Ny,Nx));
colorbar;
axis equal;
end
colormap hot;

%% 

%build projector on to lowest band
Nlowband = 4;
ProjForm = sparse(1:Nx*Ny,ones(1,Nx*Ny),psivects(:,1),Nx*Ny,Nx*Ny);
for ii=2:Nlowband
    ProjForm = ProjForm + sparse(1:Nx*Ny,ii*ones(1,Nx*Ny),psivects(:,ii),Nx*Ny,Nx*Ny);
end

%ProjForm(:,1:Nlowband) = psivects(:,1:Nlowband);
P = ProjForm*ProjForm';

X = sparse(1:Nx*Ny,1:Nx*Ny,pos(:,1));
Y = sparse(1:Nx*Ny,1:Nx*Ny,pos(:,2));
Xlowband = P*X*P;
%Ylowband = P*Y*P;
%XYcommutator = Xlowband*Ylowband-Ylowband*Xlowband;

[xeigvects,xpos] = eigs(Xlowband,Neigvects);
[xpos,sortedIndices] = sort(diag(xpos));
xeigvects = xeigvects(:,sortedIndices);
%need to convert back to x,y 2D grid

figure
for ii=1:Nlowband
subplot(3,3,ii)
imagesc(reshape(xeigvects(:,ii),Ny,Nx));
colorbar
axis equal;
end
colormap hot;
