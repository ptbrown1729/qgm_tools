function [DistanceGrid] = ellipticalGrid(ImgSize,P,UnitSystem)
%[DistanceGrid] = ellipticalGrid(ImgSize,P)
%Distance returned depends on chosen unit system. Given a point (x,y), it
%will lie on some ellipse. 'minor' unit system returns the length of the
%minor axis of that ellipse, 'major' the major axis length, and 'mean' the
%geometric mean of the two (i.e. sqrt(AB))
%ImgSize = [Ny,Nx]
%Preferred form for params: P = [Cx,Cy,AspectRatio,Theta]
%Param form allowed for backwards compatability: P = [Cx,Cy,OmegaX,OmegaY,Theta];

if ~exist('UnitSystem','var')
    UnitSystem = 'minor';
end

if ~(strcmp(UnitSystem,'minor')||strcmp(UnitSystem,'major')||strcmp(UnitSystem,'mean'))
    errorStruct.message = 'UnitSystem argument was not minor,major, or mean.';
    errorStruct.identifier = 'ellipticalGrid:argument';
end

%Check arguments...
if sum(isnan(P)) ~= 0 || sum(isinf(P)) ~= 0
    errorStruct.message = 'P contained nans or infs.';
    errorStruct.identifier = 'ellipticalGrid:nanOrinf';
    error(errorStruct);
end

if length(P)==5
    %this first case is kept for backwards compatability.
    Cx = P(1); Cy = P(2);
    AspectRatio  = P(3)/P(4);
    Theta = P(5);
else
    Cx = P(1); Cy = P(2);
    AspectRatio = P(3);
    Theta = P(4);
end

Ny = ImgSize(1); Nx = ImgSize(2);
% [X,Y] = meshgrid(1:Ny,1:Nx);
X = kron(ones(Ny,1),1:Nx); %This is equivalent to meshgrid...should use meshgrid instead...
Y = kron((1:Ny)',ones(1,Nx));
DistanceGrid = sqrt(((X-Cx)*cos(Theta)-(Y-Cy)*sin(Theta)).^2+((Y-Cy)*cos(Theta)+(X-Cx)*sin(Theta)).^2*(AspectRatio)^2);

%Quick explanation of the 'unit systems' here...
%Ellipse equation is 1 = (x-cx)^2/A^2 + (y-cy)^2/B^2
%Define d_A = sqrt((x-cx)^2+(y-cy)^2*(A/B)^2) ... which is the distance we used above
%Define d_B = sqrt((x-cx)^2*(B/A)^2 + (y-cy)^2) = (B/A)*d_A
%Define d_AB = sqrt((x-cx)^2*(B/A) + (y-cy)^2*(A/B)) = sqrt(B/A)*d_A
%for a given ellipse, d_A is the distance along the A axis, d_B along the B
%axis, and d_AB along 45 deg axis. i.e. d_A(x,y) gives the length of the A
%axis of an ellipse with the given axes A and B that contains (x,y).
%d_AB is the most natural in some sense because it gives sqrt(AB) as the
%distance...

if AspectRatio<1
    %if aspect ratio<1 we are already in 'minor' units.
    if strcmp(UnitSystem,'minor')
        
    elseif strcmp(UnitSystem,'major')
        DistanceGrid = DistanceGrid/AspectRatio;
    elseif strcmp(UnitSystem,'mean')
        DistanceGrid = DistanceGrid/sqrt(AspectRatio);
    end
elseif AspectRatio>1
    %if aspect ratio>1 we are already in 'major' units
    if strcmp(UnitSystem,'minor')
        DistanceGrid = DistanceGrid/AspectRatio;
    elseif strcmp(UnitSystem,'major')
 
    elseif strcmp(UnitSystem,'mean')
        DistanceGrid = DistanceGrid/sqrt(AspectRatio);
    end
end

end

