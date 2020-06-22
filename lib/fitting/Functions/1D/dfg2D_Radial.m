function [V,PNames,ArgStr] = dfg2D_Radial(P,X)
%[V,PNames,ArgStr] = dfg2D_Radial(P,X)
%P = [Cx,BetaMu,BetaV,Bg,Omega,TOF]
%X in SI units (meters)
%Distribution of a non-interacting 2D Fermi gas at finite temperature and
%finite time of flight.
%Best way to use this function is to fix Omega and TOF from measurements.
%BetaV = Beta*0.5*m*Omega^2

PNames = {'Cx','BetaMu','BetaV','Bg','Omega','TOF'};
ArgStr = 'pi*(mLi/h)^2*Omega^2/(1+Omega^2*TOF^2)^2*(1/BetaV)*log(1+exp(BetaMu-BetaV*(X-Cx).^2))+Bg';

K=1.3806488e-23;
mLi=9.98834e-27; %Kg
h=6.62606957e-34; %Js


if isempty(P)||isempty(X)
    V = 0;
else  
        Cx = P(1); BetaMu = P(2); BetaV = P(3); Bg = P(4); Omega = P(5); TOF = P(6);
        %Omega = 2*pi*100;
        V =  pi*(mLi/h)^2*Omega^2/(1+Omega^2*TOF^2)^2*(1/BetaV)*log(1+exp(BetaMu-BetaV*(X-Cx).^2))+Bg;
%         Cx = P(1); BetaMu = P(2); BetaV = P(3); Beta = P(4); Bg = P(5);
%         V =  2*pi*mLi/(h^2*Beta)*log(1+exp(BetaMu-BetaV*(X-Cx).^2))+Bg;
        
        
end
end
