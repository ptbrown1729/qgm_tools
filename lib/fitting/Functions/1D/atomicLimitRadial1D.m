function [V,PNames,ArgStr] = atomicLimitRadial1D(P,X)
%[V,PNames,ArgStr] = atomicLimitRadial1D(P,X)
%P = [Beta,Omega,Mu,A]
%everything in units of U, which then doesn't need to be a fit parameter
%assuming U>0 i.e. is repulsive.
%Beta -> Beta*U
%Mu -> Mu/U
%Offset such that Mu = 0 at half filling.
%Omega -> sqrt(mLi/U)*Omega if X is in real units (i.e. meters)
%Omega -> sqrt(mLi/U)*Omega*a if X is in lattice sites
%If you are using an azimuthal average in 'geometric mean units', Omega =
%sqrt(Omega_x*Omega_y).

% PNames = {'Beta','U','Omega','Mu','A'};
PNames = {'Beta','Omega','AvgMu','DeltaMu','A','SignOfU'};
ArgStr = '';

% mLi = 9.988340000000000e-27; %mass of lithium-6 in kg
% kb = 1.380650400000000e-23; %boltzmann
% h = 6.626069570000000e-34; %planck

if isempty(P)||isempty(X)
    V = 0;
elseif sum(isinf(P))||sum(isnan(P))
    V = 0;
else

    Beta = P(1); Omega = P(2); AvgMu = P(3); DeltaMu = P(4); A = P(5); SignOfU = sign(P(6));
    Pot = @(r) 0.5*Omega^2*r.^2;
    Zg = @(r) (1+2*exp(Beta*(AvgMu+0.5)).*cosh(Beta*DeltaMu).*exp(-Beta*Pot(r))+exp(2*Beta*(AvgMu+0.5)).*exp(-Beta*(2*Pot(r)+SignOfU*1)));
    V = A*(2*exp(Beta*(AvgMu+0.5))*exp(-Beta*Pot(X))+ 2*exp(2*Beta*(AvgMu+0.5)).*cosh(Beta*DeltaMu).*exp(-Beta*(2*Pot(X)+SignOfU*1)))./Zg(X);
end
end
