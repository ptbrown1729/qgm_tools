function [V,PNames,ArgStr] = ellipse(P, X, Y)
%P = [Cx,Cy,Radius,Width,Theta,Aspectratio]
%OmegaX, OmegaY
%Vary width, DeltaK to get something resonable.
PNames = {'Cx', 'Cy', 'Radius', 'Width', 'Theta', 'AspectRatio'};
ArgStr = '';

if isempty(P)||isempty(X)||isempty(Y)
    V = 0;
else
    Cx = P(1); 
    Cy = P(2); 
    Radius = P(3); 
    Width = P(4); 
    Theta = P(5); 
    AspectRatio = P(6);
    EllipseDist = sqrt(((X - Cx) .* cos(Theta) - (Y  - Cy) .* sin(Theta)).^2 +...
        (AspectRatio).^2 .* ((Y - Cy) * cos(Theta) + (X - Cx) * sin(Theta)).^2);
    ZeroIfLess = heaviside(EllipseDist - (Radius - Width / 2));
    ZeroIfGreater = heaviside((Radius + Width / 2) - EllipseDist);
    V = ZeroIfGreater .* ZeroIfLess;
    V(V~=0) = 1;
    
end
end

