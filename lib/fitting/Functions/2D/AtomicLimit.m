function [NS, densityOdd, densityEven] = AtomicLimit(z, T, U, w, plotFlag)


M = 6*1.66e-27;
kB = 1.38e-23;
hbar = 1.055e-34;

lambda = 1064e-9;
a = lambda/(2^0.5);   %lattice spacing
b = 1/(kB*T);
Er = (hbar*pi/a)^2/(2*M);
U = U*Er;

n = 35;
[X,Y] = meshgrid(1:n,1:n);
X = X-0.5*ones(n,n);
Y = Y-0.5*ones(n,n);

E = 0.5*M*w*w*(a*a*(X.*X+Y.*Y));
Z = ones(n,n)+2*z*exp(-b.*E)+z*z*exp(-2*b.*E-b*U);
Narray = 2*z*exp(-b.*E)./Z + 2*z*z*exp(-2*b.*E-b*U)./Z;
ZeroOccupancy = ones(n,n)./Z;
SingleOccupancy = 2*z*exp(-b.*E)./Z;
DoubleOccupancy = z*z*exp(-2*b.*E-b*U)./Z;
S = -(ZeroOccupancy.*log(ZeroOccupancy)+SingleOccupancy.*log(SingleOccupancy)+DoubleOccupancy.*log(DoubleOccupancy));
% Ncomplete = cat(2,flip(cat(1,flip(Narray,1),Narray),2),cat(1,flip(Narray,1),Narray));
% % SingleOccupancycomplete = cat(2,flip(cat(1,flip(SingleOccupancy,1),SingleOccupancy),2),cat(1,flipdim(SingleOccupancy,1),SingleOccupancy));
% % N2D=2*pi*U/(M*a^2*w^2)
% Scomplete = cat(2,flip(cat(1,flip(S,1),S),2),cat(1,flip(S,1),S));

% close all
% imagesc(Ncomplete)
% colorbar
% axis image

% N=4*sum(sum(Ncomplete,2),1);
% SperP = sum(sum(Scomplete))/N;

N=4*sum(sum(Narray,2),1);
SperP = 4*sum(sum(S))/N;
NS = [N SperP]';

densityOdd = SingleOccupancy(1,:);
densityEven = ZeroOccupancy(1,:) + DoubleOccupancy(1,:);

if (exist('plotFlag') == 1)
    Ncomplete = cat(2,flip(cat(1,flip(Narray,1),Narray),2),cat(1,flip(Narray,1),Narray));
    SingleOccupancycomplete = ...
        cat(2,flip(cat(1,flip(SingleOccupancy,1),SingleOccupancy),2),...
        cat(1,flip(SingleOccupancy,1),SingleOccupancy));
    Scomplete = cat(2,flip(cat(1,flip(S,1),S),2),cat(1,flip(S,1),S));
    figure
    imagesc(Ncomplete)
    colorbar
    axis image
    plot(Ncomplete(n,n+1:end),'b')
    hold on
    plot(SingleOccupancycomplete(n,n+1:end),'r')
    plot(Scomplete(n,n+1:end),'m')
    legend('number density','single occupancy','entropy')
    TforTitle = round(T*kB/Er,3,'significant');
    muForTitle = round(log(z)*kB*T/Er,3,'significant');
    title({['T = ',num2str(TforTitle),' Er'];['N = ',num2str(int16(N))];...
        ['S/particle = ',num2str(SperP),' kB'];...
        ['mu = ',num2str(muForTitle),' Er']}) 
    figure
    imagesc(Ncomplete)
    axis equal
    tilefigs
end

end