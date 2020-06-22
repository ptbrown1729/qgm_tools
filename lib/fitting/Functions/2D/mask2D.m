function [V, PNames, ArgStr] = mask2D(P, X, Y)
%[V, PNames, ArgStr] = mask2D(P, X, Y)
%P = [Cx, Cy, A, Theta, Scale, Bg, Inversion]
PNames = {'Cx', 'Cy', 'A', 'Theta', 'Scale', 'Bg', 'Inversion'};
ArgStr = '';

%create mask
SizeY = 768;
SizeX = 1024;

Radius = 100;

CenterX = SizeX / 2;
CenterY = SizeY / 2;

Circle = zeros(SizeY, SizeX);

[Xs, Ys] = meshgrid(1:size(Circle, 2), 1:size(Circle, 1));
Xs = Xs - CenterX;
Ys = Ys - CenterY;
Rs = sqrt(Xs.^2 + Ys.^2);

Circle(Rs < Radius) = 1;

Mask = Circle;
Width = 10;
Mask(floor(SizeY / 2 - Width / 2) : ceil(SizeY / 2 + Width / 2), :) = ...
    ones(size(Mask(floor(SizeY / 2 - Width / 2) : ceil(SizeY / 2 + Width / 2), :)));

Line2 = zeros(size(Mask));
Line2(floor(SizeY / 2 - Width / 2) : ceil(SizeY / 2 + Width / 2), ceil(SizeX / 2): end) = ...
    ones(size(Line2(floor(SizeY / 2 - Width / 2) : ceil(SizeY / 2 + Width / 2), ceil(SizeX / 2) : end)));
Line2 = imrotate(Line2, 30, 'bilinear', 'crop');

Mask(Line2 == 1) = 1;

%create interpolating function for mask.
[XInterp, YInterp] = meshgrid(1 : size(Mask, 2), 1 : size(Mask, 1));
XInterp = XInterp - size(Mask,2) / 2;
YInterp = YInterp - size(Mask,1) / 2;
MaskFn = @(Xq, Yq) interp2(XInterp, YInterp, Mask, Xq, Yq);
%scaled,shifted,rotated,function
FitFn = @(Xq, Yq, P) P(3) * MaskFn((Xq - P(1)) * cos(P(4)) / P(5) - (Yq - P(2)) * sin(P(4)) / P(5),...
                               (Xq - P(1)) * sin(P(4)) / P(5) + (Yq - P(2)) * cos(P(4)) / P(5));

if isempty(P) || isempty(X)
    V = 0;
else
    Cx = P(1); 
    Cy = P(2); 
    A = P(3); 
    Theta = P(4); 
    Scale = P(5); 
    Bg = P(6);
    Inversion = P(7);
    V = FitFn(X, Y, [Cx, Cy, A, Theta, Scale, Bg]);
    V(isnan(V)) = 0;
end

end
