function [V,PNames,ArgStr] = latticeFGRadial_VaryingLatticeDepth1D(P,X)
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

PNames = {'Beta','Omega','Mu','Vo','tEr','A'};
ArgStr = '';

if isempty(P)||isempty(X)
    V = 0;
elseif sum(isinf(P))||sum(isnan(P))
    V = 0;
else
    %Assume a finite size lattice system with
    %system size in each direction (# of sites);
    N = 50;
    %I think the way you want to think about setting the size is:
    %you need the separation between the allowed k vectors in energy to be
    %small compared with the temperature. Roughly, BandWidth/N=8/N<<T
    
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
    
    
    
    
    Beta = P(1); Omega = P(2); Mu = P(3); Vo = P(4); tEr = P(5); A = P(6);
    %for our lattice, approximation for tunneling in terms of depth is
    %t = 1.5*Depth*exp(-2*sqrt(Depth)), in units of Er
    %and t(r) = t(r)/to
    t = @(r) (1-0.5*Omega^2*r.^2/Vo).*exp(-2*sqrt(tEr)*sqrt(Vo-0.5*Omega^2*r.^2))/exp(-2*sqrt(tEr)*sqrt(Vo));
    
    Es = -2*(cos(Kxxs)+cos(Kyys)).*t(XX);
    Mus = (Mu - 0.5*Omega.^2*XX.^2); %.*t(XX);
%     Betas = Beta./t(XX);
    Occ = 1./(exp(Beta.*(Es - Mus))+1);
    %now compute density/site by summing over allowed states
    V = A*sum(sum(Occ,ndims(XX)-1),ndims(XX))/N^2;
    V(abs(imag(V))>1e-7) = 0;
    V(abs(imag(V))>1e-3*abs(real(V))) = 0;
    V = real(V);
    
end

end