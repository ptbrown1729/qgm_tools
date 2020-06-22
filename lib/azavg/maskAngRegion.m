function [MaskedImg] = maskAngRegion(P,Img)
%maskAngRegion
Cx = P(1); Cy = P(2); Theta1 = P(3); Theta2 = P(4);

if length(P) == 5
    R1 = P(5); 
    %R2 = P(6);
else
    R1  = 0; R2 = sqrt(size(Img,1)^2+size(Img,2)^2)+1;
end

%if Cx or Cy is an integer, will cause problems in deciding if a pixel is
%to the right or left of the center.
if mod(Cx,1) == 0
    Cx = Cx + 1e-5;
end
if mod(Cy,1) == 0
    Cy = Cy +1e-5;
end

[Ny,Nx] = size(Img);
Xpts = kron(ones(Ny,1),1:Nx);
Ypts = kron((1:Ny)',ones(1,Nx));
DistanceArray = sqrt((Xpts-Cx).^2+(Ypts-Cy).^2);
AngleArray = atan(-(Ypts-Cy)./(Xpts-Cx));

%Can't distinguish theta and theta+pi with atan. Need additional matrix to
%distinguish this. Zero if less than pi, one if more.
ThetaPlusPi = ones(Ny,Nx);
ThetaPlusPi(Cx<=Xpts)=0;

Theta1 = mod(Theta1,2*pi);
Theta2 = mod(Theta2,2*pi);
if Theta1>Theta2
    ThetaLess = Theta2;
    Theta2 = Theta1;
    Theta1 = ThetaLess;
end

%get Thetas in range [-pi/2,pi/2] and proper quadrant
if (0 <= Theta1) && (Theta1 < pi/2)
    HalfBool1 = 0;
elseif (3*pi/2 < Theta1) && (Theta1 <= 2*pi)
    HalfBool1 = 0;
    Theta1 = Theta1-2*pi;
else
    HalfBool1 = 1;
    Theta1 = Theta1 - pi;
end

if (0 <= Theta2 ) && (Theta2 < pi/2)
    HalfBool2 = 0;
elseif (3*pi/2 < Theta2) && (Theta2 <= 2*pi)
    HalfBool2 = 0;
    Theta2 = Theta2-2*pi;
else
    HalfBool2 = 1;
    Theta2 = Theta2 - pi;
end

%Create mask
if HalfBool1 == HalfBool2 && Theta1<Theta2
    Mask = ones(Ny,Nx);
    Mask((Theta1<AngleArray)&(AngleArray<Theta2)) = 0;
    Mask(ThetaPlusPi ~= HalfBool1) = 1;
elseif HalfBool1 == HalfBool2 && Theta2<Theta1 %wraps around.
    Mask = ones(Ny,Nx);
    Mask(Theta1<AngleArray) = 0;
    Mask(Theta2>AngleArray) = 0;
    Mask(ThetaPlusPi ~= HalfBool1) = 0;
else
    Mask1 = ones(Ny,Nx);
    Mask1((Theta1<AngleArray)&(AngleArray<=pi/2+1e-5)) = 0;
    %Mask1(ThetaPlusPi = 1) = 1;
    Mask1(ThetaPlusPi ~= HalfBool1) = 1;
    Mask2 = ones(Ny,Nx);
    Mask2((-pi/2<AngleArray)&(AngleArray<Theta2)) = 0;
    Mask2(ThetaPlusPi ~= HalfBool2) = 1;
    Mask = Mask1.*Mask2;    
end

Mask(DistanceArray<R1) = 1;
%Mask(DistanceArray>R2) = 1;

MaskedImg = Img.*Mask;




end

