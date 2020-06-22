function [V, PNames, ArgStr] = gaussian2D(P, X, Y)
% Rotated 2D Gaussian
%
% [V] = [V,PNames,ArgStr] = gaussian2D(P,X,Y)
% P = [CxG, CyG, SxG, SyG, AmpG, ThetaG, Bg]
% Theta is the counter-clockwise angle between the axis with width SxG, if
% we regard the y axis as pointing downards (i.e. as matlab would display a
% matrix).

PNames = {'CxG', 'CyG', 'SxG', 'SyG', 'AG', 'ThetaG', 'Bg'};
ArgStr = 'Bg+AG*(exp(-0.5*((X-CxG)*cos(ThetaG)-(Y-CyG)*sin(ThetaG)).^2./(SxG^2)-0.5*((Y-CyG)*cos(ThetaG)+(X-CxG)*sin(ThetaG)).^2./(SyG^2)))';

if isempty(P)||isempty(X)||isempty(Y)
    V = 0;
else
    CxG = P(1); 
    CyG = P(2); 
    SxG = P(3); 
    SyG = P(4); 
    AG = P(5); 
    ThetaG = P(6); 
    Bg = P(7);
    V = Bg + AG * (exp(-0.5 * ( (X - CxG) * cos(ThetaG) - (Y - CyG) * sin(ThetaG)).^2 ./ (SxG^2) -...
        0.5 * ( (Y-CyG) * cos(ThetaG) + (X - CxG) * sin(ThetaG)).^2 ./ (SyG^2)));
end
end

