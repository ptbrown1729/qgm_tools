function [V,PNames,ArgStr] = Fermi1DModConf(P,X)
%P = [Temp, Depth, mu, beamRadius, bg]
%NOTE: None of the arguments can be single ascii letters.
PNames = {'T','Depth','mu','beamRadius','bg'};
ArgStr = '2*pi*9.98834e-27*1.3806488e-23*T/(6.62606957e-34^2)*log(1+exp((Depth*exp(-2*x.*x/(beamRadius*beamRadius))+mu)/(1.3806488e-23*T)))+bg';

kB=1.3806488e-23;
mLi=9.98834e-27; %Kg
h=6.62606957e-34; %Js
pixelSize=8.667e-7; %m

if isempty(P)||isempty(X)
    V = 0;
else
    T = P(1); Depth = kB*P(2); mu = kB*P(3); beamRadius = P(4); bg=P(5); X = X*pixelSize;
    V = 2*pi*mLi*kB*T/(h^2)*log(1+exp((Depth*(exp(-2*X.*X/(beamRadius*beamRadius))-1)+mu)/(kB*T)))+bg;
end

end

