function [V,PNames,ArgStr] = latticeFGRadial_GaussPot1D(P,X)
%[V,PNames,ArgStr] = latticeFGRadial1D(P,X)
%Non-interacting FG in lattice, in tight-binding limit and LDA.
%All parameters in units of t, the tunneling
%chemical potential is referenced so that mu=0 is half filling.
%Dispersion is E(k) = -2*(cos(kxa)+cos(kya)).
%P = [Beta,Waist,Potential0,Mu,A]
%Beta = Beta*t
%Mu = Mu/t
%Waist is in the same units as X.
%Omega = sqrt(2*Pot0/Waist^2) in the same units...
%To convert Omega to real units...
%Omega = sqrt(mLi/t)*Omega if X is in real units
%Omega = sqrt(mLi/t)*a*Omega if X is in 'lattice sites'. Where a is the
%lattice spacing.

%If you are using an azimuthal average in 'geometric mean units', Omega =
%sqrt(Omega_x*Omega_y).

PNames = {'Beta','Waist','Pot0','Mu','A'};
ArgStr = '';

if isempty(P)||isempty(X)
    V = 0;
elseif sum(isinf(P))||sum(isnan(P))
    V = 0;
else
    %Assume a finite size lattice system with
    %system size in each direction (# of sites);
    N = 30;
    %period boundary conditions imply the allowed momenta are
    kx = (2*pi)/(N)*(0:N-1);
    ky = (2*pi)/(N)*(0:N-1);
    %put them in a grid...
    [Kxs,Kys] = meshgrid(kx,ky);
    %These are all of the allowed states, and we know there energies, so the
    %density is just going to be
    %\int f(E(k))*A/(2*pi)^2 dkxdky -> sum_{n,m} f(-2tcos(2pi*m/N)-2tcos(2pi*n/N))
    %so when we convert to the sum, we lose the dependence on system size...
    
    %%%need to vectorize to accept arbitrary X shapes...
    %do this by repeating X over the last two dimensions
    %and repeating Kxs and Kys over the first n-dimensions.
    %later will sum over Kxs and Kys by summing the last two dimensions.
    XX = repmat(X,[ones(1,ndims(X)),length(ky),length(kx)]);
    Kxxs = permute(repmat(Kxs,[1,1,size(X)]),[3:ndims(X)+2,1,2]);
    Kyys = permute(repmat(Kys,[1,1,size(X)]),[3:ndims(X)+2,1,2]);
    
    Beta = P(1); Waist = P(2); Pot0 = P(3); Mu = P(4); A = P(5);
    Es = -2*(cos(Kxxs)+cos(Kyys));
    Mus = Mu - Pot0*(1-exp(-2*XX.^2/Waist^2));
    Occ = 1./(exp(Beta*(Es - Mus))+1);
    %now compute density/site by summing over allowed states
    V = A*sum(sum(Occ,ndims(XX)-1),ndims(XX))/N^2;
    
end

end