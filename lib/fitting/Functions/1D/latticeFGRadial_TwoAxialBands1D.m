function [V,PNames,ArgStr] = latticeFGRadial_TwoAxialBands1D(P,X)
%[V,PNames,ArgStr] = latticeFGRadial1D(P,X)
%Non-interacting FG in lattice, in tight-binding limit and LDA.
%All parameters in units of t, the tunneling
%chemical potential is referenced so that mu=0 is half filling.
%Dispersion is E(k) = -2*(cos(kxa)+cos(kya)).
%P = [Beta,Omega,Mu,A]
%Beta = Beta*t
%Mu = Mu/t
%Omega = sqrt(mLi/t)*Omega if X is in real units
%Omega = sqrt(mLi/t)*a*Omega if X is in 'lattice sites'. Where a is the
%lattice spacing.
%If you are using an azimuthal average in 'geometric mean units', Omega =
%sqrt(Omega_x*Omega_y).

PNames = {'Beta','Omega','Mu','BandGap','A'};
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
    
    Beta = P(1); Omega = P(2); Mu = P(3); BandGap = P(4);
    A = P(5);
    Es = -2*(cos(Kxxs)+cos(Kyys));
    Mus = Mu - 0.5*Omega.^2*XX.^2;
    OccLowerBand = 1./(exp(Beta*(Es - Mus))+1);
    %extra factor of 4, i.e. half the bandwidth, so that BandGap is really
    %the BandGap and not just some energy offset.
    %Note that this is the second AXIAL band, and not the second lattice
    %band, which would have a different dispersion...
    OccUpperBand = 1./(exp(Beta*(Es + (BandGap+4) - Mus))+1);
    %now compute density/site by summing over allowed states
    %The minus sign here is the parity imaging coming into play...
    %Assuming that T<<BandGap, the upper band only begins to matter after
    %the lower band is completely filled. Then, we will always see an atom
    %on a site, unless the higher band is occupied. So the average filling
    %from the higher band should be subtracted from the lower band...
    
    V = A*sum(sum(OccLowerBand-OccUpperBand,ndims(XX)-1),ndims(XX))/N^2;
    
end

end