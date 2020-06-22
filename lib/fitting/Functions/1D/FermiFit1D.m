function [V,PNames,ArgStr] = FermiFit1D(P,X)
%[V,PNames,ArgStr] = FermiFit1D(P,X)
%Natural energy scale is the Fermi-Energy, i.e. mu
%P = [BetaMu,Omega,Beta,A]
%BetaMu = Beta*Mu, i.e. Beta in units of Mu
%Omega = sqrt(mLi/Mu)*Omega
%Beta

%P = [BetaMu,BetaV Omega]
% PNames = {'CG','SG','AG'};
% ArgStr = '';

K=1.3806488e-23;
mLi=9.98834e-27; %Kg
h=6.62606957e-34; %Js

if isempty(P)||isempty(X)
    V = 0;
else
    BetaMu = P(1); Omega = P(2); Beta = P(3);
    V = A*2*pi*mLi/(h^2*Beta)*log(1+exp(BetaMu*(1-Omega^2*X.^2)));
    
%     BetaMu = P(1); BetaV = P(2); Beta = P(3);
%     V = 2*pi*mLi/(h^2*Beta)*log(1+exp(BetaMu-BetaV*X.^2));
%     V = pi*(mLi*Omega/h)^2*(1/BetaV)*log(1+exp(BetaMu-BetaV*X.^2));
end

end

