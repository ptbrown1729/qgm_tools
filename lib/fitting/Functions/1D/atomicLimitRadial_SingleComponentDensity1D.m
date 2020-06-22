function [V,PNames,ArgStr] = atomicLimitRadial_SingleComponentDensity1D(P,X)
%[V,PNames,ArgStr] = atomicLimitRadial_Doubles1D(P,X)
%X can be in real distance units or lattice sites.
%P = [Beta,Omega,AvgMu,DeltaMu,A,SignOfU]
%everything in units of U, which then doesn't need to be a fit parameter
%Beta -> Beta*U
%AvgMu -> AvgMu/U
%DeltaMu -> DeltaMu/U
%Chose offset so that half-filling occurs at AvgMu = 0
%Omega -> sqrt(mLi/U)*Omega if X is in real units
%Omega -> sqrt(mLi/U)*a*Omega if X is in 'lattice sites'
%Where a is the lattice spacing.
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
    %rewrite partition function in terms of AvgMu=0.5*(mu_1+mu_2) and
    %DeltaMu=0.5*(mu_1-mu_2)...
    Zg = @(r) (1+2*exp(Beta*(AvgMu+0.5)).*cosh(Beta*DeltaMu).*exp(-Beta*Pot(r))+exp(2*Beta*(AvgMu+0.5)).*exp(-Beta*(2*Pot(r)+SignOfU*1)));
%     V = A*(2*exp(Beta*AvgMu).*cosh(Beta*DeltaMu).*exp(-Beta*Pot(X)))./Zg(X);
    V = A*(exp(Beta*(AvgMu+0.5)).*cosh(Beta*DeltaMu).*exp(-Beta*Pot(X))+exp(2*Beta*(AvgMu+0.5)).*exp(-Beta*(2*Pot(X)+SignOfU*1)))./Zg(X);
end
end
