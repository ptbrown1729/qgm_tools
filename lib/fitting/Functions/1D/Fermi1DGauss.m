function [V,PNames,ArgStr] = Fermi1DGauss(P,X)
%P = [Temp, Depth, mu, beamRadius, bg]
%NOTE: None of the arguments can be single ascii letters.
PNames = {'T','Omega','mu','Sigma','bg'};
ArgStr = '2*pi*9.98834e-27*1.3806488e-23*T/(6.62606957e-34^2)*log(1+exp((Depth*exp(-2*x.*x/(beamRadius*beamRadius))+mu)/(1.3806488e-23*T)))+bg';

kB=1.3806488e-23;
mLi=9.98834e-27; %Kg
h=6.62606957e-34; %Js
% pixelSize=8.667e-7; %m
Density2OD = 3*671e-9^2/(2*pi);

if isempty(P)||isempty(X)
    V = 0;
else
    T = P(1); Omega = P(2); mu = kB*P(3); Sigma = P(4); bg=P(5);
%     0.25*mLi*Omega^2*Sigma^2
    V = Density2OD*2*pi*mLi*kB*T/(h^2)*log(1+exp((0.25*mLi*Omega^2*Sigma^2*(exp(-2*X.*X/Sigma^2)-1)+mu)/(kB*T)))+bg;
end

end