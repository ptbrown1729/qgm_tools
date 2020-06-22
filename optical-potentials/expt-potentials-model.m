%constants
K=1.3806488e-23;
mLi=9.98834e-27; %Kg
h=6.62606957e-34; %Js
hbar=h/(2*pi); %Js
c=2.99792458e8; %m/s
e=1.602176565e-19; %C
me=9.10938291e-31; %kg
mp=1.672621777e-27; %kg

%6Li D2 Line
Lambda0=670.977e-9; %nm
Omega0 =2*pi*c/Lambda0;
Gamma=2*pi*5.8724e6; %Hz
Isat=25.4; %W/m^2
%% Light Sheet

LambdaLS = 1070e-9;
OmegaLS = 2*pi*c/LambdaLS;

P_LS = 300e-3;
W0Z_LS = 5.1e-6;
W0X_LS = 49e-6;
%Derived
ZX_LS = pi*W0X_LS^2/LambdaLS;
WX_LS = @(Y) W0X_LS*sqrt(1+(Y/(ZX_LS)).^2);
ZZ_LS = pi*W0Z_LS^2/LambdaLS;
WZ_LS = @(Y) W0Z_LS*sqrt(1+(Y/(ZZ_LS)).^2);
%Potential
Alpha_LS = -3*pi*c^2/(2*Omega0^3)*(Gamma/(Omega0 - OmegaLS) + Gamma/(Omega0 + OmegaLS));
V_LS = @(X,Y,Z)Alpha_LS*(2*P_LS./(pi*WX_LS(Y).*WZ_LS(Y))).*exp(-2*(X.^2./WX_LS(Y).^2+Z.^2./WZ_LS(Y).^2));

%Trapping Frequencies and Depth
Depth_Atoms = abs(V_LS(0,0,0)/K);
FX0 = sqrt(4*abs(V_LS(0,0,0))/(mLi*W0X_LS^2))/(2*pi);
FZ0 = sqrt(4*abs(V_LS(0,0,0))/(mLi*W0Z_LS^2))/(2*pi);
FY0 = sqrt(2*abs(V_LS(0,0,0))*0.5*(1/ZZ_LS^2+1/ZX_LS^2)/mLi)/(2*pi);

V_Harm_LS = @(X,Y,Z) 0.5*mLi*((2*pi*FX0)^2*X.^2+(2*pi*FY0)^2*Y.^2+(2*pi*FZ0)^2*Z.^2)+V_LS(0,0,0);

VDevX = @(X) 0.1-abs(V_LS(X,0,0)-V_Harm_LS(X,0,0))/abs(V_LS(X,0,0));
%RX_HarmoDev = fnzeros(VDevX,[0,2*W0X_LS])

disp(cat(2,num2cell([Depth_Atoms,FX0,FZ0,FY0])',{'Depth Atoms','FX0 LS','FY0 LS','FZ0 LS'}'))

X = linspace(-W0X_LS*2,W0X_LS*2,100);
plot(X,V_Harm_LS(X,0,0),'r')
hold on;
plot(X,V_LS(X,0,0),'b.')
hold off;

%% Pancakes

LambdaPC = 532e-9;
KPC = 2*pi/LambdaPC;
OmegaPC = c*KPC;

P_PC = 0.2; %for a single beam.
W0Z_PC = 25e-6;
W0X_PC = 120e-6;
ThV = 7*pi/180;
ThH = 135*pi/180;

%Derived
ZX_PC = pi*W0X_PC^2/LambdaPC;
WX_PC = @(Y) W0X_PC*sqrt(1+(Y/(ZX_PC)).^2);
ZZ_PC = pi*W0Z_PC^2/LambdaPC;
WZ_PC = @(Y) W0Z_PC*sqrt(1+(Y/(ZZ_PC)).^2);
%
Alpha_PC = -3*pi*c^2/(2*Omega0^3)*(Gamma/(Omega0 - OmegaPC) + Gamma/(Omega0 + OmegaPC));
%Potential. Assuming beam travels along Y in its own coordinates.
E_PC_Beam = @(X,Y,Z) sqrt((2*P_PC./(pi.*WX_PC(Y).*WZ_PC(Y))).*exp(-2*(X.^2./WX_PC(Y).^2+Z.^2./WZ_PC(Y).^2))).*exp(1i*KPC*Y);
%Euler angles to rotate beams to correct orientation.
%Transformation matrix, for reference.
Uinv = @(Theta,Phi)[cos(Phi)  sin(Phi) 0; -sin(Phi) cos(Phi) 0; 0 0 1]*[1 0 0; 0 cos(Theta) sin(Theta); 0 -sin(Theta) cos(Theta)];
E_PC_Upper = @(X,Y,Z) E_PC_Beam(X*cos(ThH)+Y*sin(ThH)*cos(ThV)+Z*sin(ThH)*sin(ThV),-X*sin(ThH)+Y*cos(ThH)*cos(ThV)+Z*cos(ThH)*sin(ThV),-Y*sin(ThV)+Z*cos(ThH));
E_PC_Lower = @(X,Y,Z)  E_PC_Beam(X*cos(ThH)+Y*sin(ThH)*cos(-ThV)+Z*sin(ThH)*sin(-ThV),-X*sin(ThH)+Y*cos(ThH)*cos(-ThV)+Z*cos(ThH)*sin(-ThV),-Y*sin(-ThV)+Z*cos(ThH));
V_PC = @(X,Y,Z) Alpha_PC*real((E_PC_Upper(X,Y,Z)+E_PC_Lower(X,Y,Z)).*conj(E_PC_Upper(X,Y,Z)+E_PC_Lower(X,Y,Z)));

%Something wrong with frequency calcs...
FZ_Central = sqrt(0.1*Alpha_PC*abs(E_PC_Beam(0,0,0))^2*KPC^2*sin(ThV)^2/mLi);
FY_OwnCoords_Single = -sqrt(2*abs(-Alpha_PC*abs(E_PC_Upper(0,0,0))^2)*0.5*(1/ZZ_PC^2+1/ZX_PC^2)/mLi)/(2*pi);
FX_OwnCoords_Single = -sqrt(4*abs(-Alpha_PC*abs(E_PC_Upper(0,0,0))^2/(mLi*W0X_PC^2)))/(2*pi);

disp(cat(2,num2cell([FX_OwnCoords_Single,FY_OwnCoords_Single,FZ_Central])',{'FX0 PC','FY0 PC','FZ0 PC'}'))

V_Harm_PC = @(X,Y,Z) -0.5*mLi*(2*pi*FZ_Central)^2*Z.^2+V_PC(0,0,0);

Z = linspace(-0.5e-6,0.5e-6,100);
plot(Z,V_PC(0,0,Z),'b.')
hold on;
plot(Z,V_Harm_PC(0,0,Z),'r');
hold off;

%% Total Potential

V_Tot = @(X,Y,Z) V_PC(X,Y,Z)+V_LS(X,Y,Z);

Y = linspace(-300e-6,300e-6,4500);
Z = linspace(-100e-6,100e-6,1500);
[y,z] = meshgrid(Y,Z);
subplot(2,2,1)
imagesc(V_Tot(0,y,z)/K*1e6)
title('Y-Z Plane')
colorbar

X = linspace(-300e-6,300e-6,4500);
Z = linspace(-100e-6,100e-6,1500);
[x,z] = meshgrid(X,Z);
subplot(2,2,2)
imagesc(V_Tot(x,0,z)/K*1e6)
title('X-Z Plane')
colorbar

X = linspace(-300e-6,300e-6,4500);
Y = linspace(-300e-6,300e-6,4500);
[x,y] = meshgrid(X,Y);
subplot(2,2,3)
imagesc(V_Tot(x,y,0)/K*1e6)
title('X-Y Plane')
colorbar

figure
subplot(2,2,1)
Z = linspace(-20e-6,20e-6,500);
plot(Z,V_Tot(0,0,Z)/K*1e6)
title('Line Along (0,0,Z)')

subplot(2,2,2)
X = linspace(-100e-6,100e-6,500);
plot(X,V_Tot(X,0,0)/K*1e6)
title('Line Along (X,0,0)')

subplot(2,2,3)
Y = linspace(-100e-6,100e-6,500);
plot(Y,V_Tot(0,Y,0)/K*1e6)
title('Line Along (0,Y,0)')

