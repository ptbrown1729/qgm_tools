function [V,PNames,ArgStr] = dfg2D(P,X,Y)
% 2D degenerate ideal Fermi gas in Free space, including time-of-flight
%
% Some of these parameters should be fixed before attempting to do any
% fitting. See, for example, the documentation on fit2D for how to fix
% these.
% P = [Cx, Cy, BetaMu, BetaVx, BetaVy, OmegaX, omegaY, TOF, Theta, Bg]
%
% X and Y must be in SI units to get Temp in K (i.e. not in 'pixels' but in
% meters)

PNames = {'Cx','Cy','BetaMu','BetaVx','BetaVy','OmegaX','OmegaY','TOF','Theta','Bg'};
ArgStr = '';

K=1.3806488e-23;
mLi=9.98834e-27; %Kg
h=6.62606957e-34; %Js

if isempty(P)||isempty(X)||isempty(Y)
    V = 0;
else
    Cx = P(1); 
    Cy = P(2); 
    BetaMu = P(3); 
    BetaVx = P(4); 
    BetaVy = P(5); 
    OmegaX = P(6); 
    OmegaY = P(7);
    TOF = P(8); 
    Theta = P(9); 
    Bg = P(10);
 
    V =  pi * (mLi / h)^2 * OmegaX * OmegaY / (1 + OmegaX^2 * TOF^2) / (1 + OmegaY^2 * TOF^2) *...
        (1 / sqrt(BetaVx * BetaVy)) *...
        log(1 + exp(BetaMu - BetaVx*(X - Cx).^2 - BetaVy * (Y - Cy).^2)) + Bg;
    
end
end

