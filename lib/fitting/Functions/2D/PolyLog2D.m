function [V,PNames,ArgStr] = PolyLog2D(P,X,Y)
%P = [Cx,Cy,Sx,Sy,Amp,Theta,Bg,Fugacity]
%Theta is measured CCW.
%Rotated 2D PolyLog to get T/TF of a 3D oblate cloud
%Expression in Pg 260 of DeMarco's thesis
%(https://jila.colorado.edu/jin/sites/default/files/files/2001_demarco.pdf)
% run
PNames = {'Cx','Cy','Sx','Sy','A','Theta','Bg','Fugacity'};
ArgStr = '';

if isempty(P) || isempty(X) || isempty(Y)
    V = 0;
else
    Cx = P(1); 
    Cy = P(2); 
    Sx = P(3); 
    Sy = P(4); 
    A = P(5); 
    Theta = P(6); 
    Bg = P(7); 
    F = P(8);
    V = Bg + A * PolyLog2(F * (exp(-0.5 * ((X - Cx) * cos(Theta) - (Y - Cy) * sin(Theta)).^2 ./ (Sx^2) -...
                                    0.5 * ((Y - Cy) * cos(Theta) + (X - Cx) * sin(Theta)).^2 ./ (Sy^2))))...
                                    /PolyLog2(F);
end
end