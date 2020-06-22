% well = @(x,y) 0.5*(x.^2+y.^2-1).*(abs(x.^2+y.^2)<1);
% V=@(x,y)0.5*(x.^2+y.^2);
V=@(x,y) 0.5*(x.^2+y.^2);
%create x-grid
xstart = -10;
xend = 10;
dx = 1/10;
xpts = xstart:dx:xend;
Nx = length(xpts);
%create y-grid
ystart = -10;
yend = 10;
dy = 1/10;
ypts = ystart:dy:yend;
Ny = length(ypts);
%build position vector
pos = [kron(xpts',ones(length(ypts),1)),kron(ones(length(xpts),1),ypts')];
%(x1,y1),(x1,y2),...(x1,yn),(x2,y1),...
%Build Hamiltonian

Vmatr = sparse(1:Nx*Ny,1:Nx*Ny,V(pos(:,1),pos(:,2)));

DxForm = -2*sparse(1:Nx,1:Nx,ones(1,Nx))+sparse(1:Nx-1,2:Nx,ones(1,Nx-1),Nx,Nx)+sparse(2:Nx,1:Nx-1,ones(1,Nx-1),Nx,Nx);
DxSqrd = kron(DxForm,sparse(1:Ny,1:Ny,ones(1,Ny)))/dx^2;
DyForm = -2*sparse(1:Ny,1:Ny,ones(1,Ny))+sparse(1:Ny-1,2:Ny,ones(1,Ny-1),Ny,Ny)+sparse(2:Ny,1:Ny-1,ones(1,Ny-1),Ny,Ny);
DySqrd= kron(sparse(1:Nx,1:Nx,ones(1,Nx)),DyForm)/dy^2;

%DxForm = sparse(diag(2*ones(Ny,1))+diag(ones(Ny-1,1),1)+diag(ones(Ny-1,1),-1));
%DxSqrd = sparse(kron(DxForm,diag(ones(Nx,1))))/dx^2;
%DyForm = sparse(diag(2*ones(Nx,1))+diag(ones(Nx-1,1),1)+diag(ones(Nx-1,1),-1));
%DySqrd = sparse(kron(diag(ones(Ny,1)),DyForm))/dy^2;

H = -0.5*(DxSqrd+DySqrd) + Vmatr; %(hbar^2/(2*mLi))*
%diagonalize Hamiltonian
Neigvects = 7;
[psivects,E] = eigs(H,Neigvects,'lm'); %note 'sm' is "smallest magnitude"
[E,sortedIndices] = sort(diag(E));
psivects = psivects(:,sortedIndices);

subplot(3,4,1)
imagesc(reshape(V(pos(:,1),pos(:,2)),Ny,Nx));
colorbar;
axis equal;

subplot(3,4,2)
plot([0,1],[E,E],'b')

for ii=1:7
subplot(3,4,ii+2)
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
