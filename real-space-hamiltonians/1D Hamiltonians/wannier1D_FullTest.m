%Full test of Wannier code. Some of this code is implemented in the function
%wannier1D.

%diagonalize Hamiltonian
%define potential
lambda = 1064e-9;
k=2*pi/lambda;
Er=hbar^2*k^2/(2*mLi);
V=@(x) 0.2*Er*cos(2*k*x);
%repeats with period lambda/2

%initialize grid
xstart = 0;
xend = 51*lambda/2-lambda/100;
dx = lambda/100;
Nx = round((xend-xstart)/dx+1);
xpts = xstart:dx:xend;

%define Hamiltonian
Potential = sparse(1:Nx,1:Nx,V(xpts));
%Potential = diag(V(xpts));
%Second derivative operator
dxSqrPeriodic =(1/dx^2)*(-2*sparse(1:Nx,1:Nx,ones(1,Nx))+sparse(1:Nx-1,2:Nx,ones(1,Nx-1),Nx,Nx)+sparse(2:Nx,1:Nx-1,ones(1,Nx-1),Nx,Nx));
H = -hbar^2/(2*mLi)*dxSqrPeriodic + Potential;

[psivects,Es]=getDiagH(Hlattice,60);


%build projector on to lowest band
Nlowband = 51;
ProjForm = sparse(1:Nx,ones(1,Nx),psivects(:,1),Nx,Nx);
for ii=2:Nlowband
    ProjForm = ProjForm + sparse(1:Nx,ii*ones(1,Nx),psivects(:,ii),Nx,Nx);
end
%ProjForm(:,1:Nlowband) = psivects(:,1:Nlowband);
P = ProjForm*ProjForm';
X = sparse(1:Nx,1:Nx,xpts);
Xlowband = P*X*P;
[wanns,xpos] = eigs(Xlowband,Nlowband);
%create projector

% m = 51;
% ProjForm = zeros(Nx,Nx);
% ProjForm(:,1:m) = psivects(:,1:m);
% Proj = ProjForm*ProjForm';
% X=diag(xpts);
% %diagonalize projector and recover wannier states
% [wann,pos] = eig(Proj*X*Proj);
% pos = diag(pos);
% pos = pos(1:m);

% %periodic case. x(n+1)=x1
% dxSqrPeriodic = dxSqr;
% dxSqrPeriodic(1,N) = 1/xspace^2;
% dxSqrPeriodic(N,1) = 1/xspace^2;
% Hperiodic = -0.5*dxSqrPeriodic + Potential;
% [psivectsPeriodic,EPeriodic] = eig(Hperiodic);
% [EPeriodic,position] = sort(diag(EPeriodic));
% mm = 6;
% ProjFormPeriodic = zeros(N,N);
% ProjFormPeriodic(:,1:mm) = psivectsPeriodic(:,1:mm);
% ProjPeriodic = ProjFormPeriodic*ProjFormPeriodic';
% XPeriodic = diag(exp(2*pi*1i*xpts/(xend-xstart+xspace)));
% [wannPeriodic,posPeriodic] = eig(ProjPeriodic*XPeriodic*ProjPeriodic);
% posPeriodic = diag(posPeriodic);
% posPeriodic = posPeriodic(1:mm);


%plot
figure
subplot(2,2,1);
plot(xpts,V(xpts))
ylabel('Potential')

subplot(2,2,2);
plot([1,2],[Es(1:60),Es(1:60)],'b')
ylabel('Energy Levels')

subplot(2,2,3);
plot(xpts,psivects(:,1:3));
ylabel('Eigenstates')

subplot(2,2,4);
plot(xpts,wanns(:,1:51))
ylabel('Wannier States')


