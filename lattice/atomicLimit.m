%atomic limit calculation.
%assuming U>>t. Treat each site as indepenent, with same temperature and
%chemical potential from LDA

%fundamental constants
h = 6.6261e-34;
k = 1.3807e-23;
mLi = 9.9883e-27;

%temperature and chemical potential
T = 90e-9;
Beta = 1/(k*T);
mu = h*6e3;

%experimntal parameters
klatt = 2*pi/(1064-9*sqrt(2));
a = 532e-9*sqrt(2);
Waist = 60e-6;
U = h*10e3;
Omega = 2*pi*500;

%Radial calculation
V = @(r) 0.5*mLi*Omega^2*r.^2;
Zg = @(r) (1+2*exp(Beta*mu)*exp(-Beta*V(r))+exp(2*Beta*mu)*exp(-Beta*(-U+2*V(r))));
NSite = @(r) (exp(Beta*mu)*exp(-Beta*V(r))+ 2*exp(2*Beta*mu)*exp(-Beta*(-U+2*V(r))))./Zg(r);
n1D = @(r) NSite(r)/a^2;

NSiteSingleSpecies = @(r)(exp(Beta*mu)*exp(-Beta*V(r))+ exp(2*Beta*mu)*exp(-Beta*(-U+2*V(r))))./Zg(r);
nSingleSpecies = @(r) NSiteSingleSpecies(r)/a^2;

Rs = -3*Waist:a:3*Waist;
Rpos = Rs(Rs>=0);
Nfrom1D = sum(2*pi*NSite(Rpos)/a^2.*Rpos*a);
%fprintf('Total atom number N = %0.0f, mu = %0.2f KHz, T = %0.2f nK \n',Nfrom1D,mu/h/1e3,T/1e-9);
hfig = figure;
subplot(3,1,1)
plot(Rs,NSite(Rs));
ylim([0,2.2])
hold on;

subplot(3,1,2)
plot(Rs,nSingleSpecies(Rs));

%2D calculation
OmegaX = 2*pi*500;
OmegaY = 2*pi*500;
V2D = @(X,Y) 0.5*mLi*(OmegaX^2*X.^2+OmegaY^2*Y.^2);
Zg2D = @(X,Y)(1+2*exp(Beta*mu)*exp(-Beta*V2D(X,Y))+exp(2*Beta*mu)*exp(-Beta*(-U+2*V2D(X,Y))));
NSite2D = @(X,Y)(exp(Beta*mu)*exp(-Beta*V2D(X,Y))+ 2*exp(2*Beta*mu)*exp(-Beta*(-U+2*V2D(X,Y))))./Zg2D(X,Y);
n2d = @(X,Y) NSite2D(X,Y)/(a^2);

Xs = -3*Waist:a:3*Waist;
Ys = -3*Waist:a:3*Waist;
[X,Y] = meshgrid(Xs,Ys);

N = sum(sum(NSite2D(X,Y)));
fprintf('Total atom number N = %0.0f, mu = %0.2f KHz, T = %0.2f nK \n',N,mu/h/1e3,T/1e-9);
subplot(3,1,3)
imagesc(NSite2D(X,Y));
axis equal
axis image
colorbar

name = sprintf('N = %0.0f, mu = %0.2f KHz, T = %0.2f nK \n',N,mu/h/1e3,T/1e-9);
set(hfig,'name',name);
%hfig.name = name;