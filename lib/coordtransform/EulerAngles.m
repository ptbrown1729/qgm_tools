function [U] = EulerAngles(Phi, Theta, Psi)
%
% %%Arguments%%
%
% Phi
%
% Theta
%
% Psi
%
% %%Returns%%
% 
% U is a 3x3 matrix which acts on a vector in R^3.
% See affine3d also!

Uz = [1 0 0; 0 cos(Phi) -sin(Phi); 0 sin(Phi) cos(Phi)];
Ux = [cos(Theta)  -sin(Theta) 0; sin(Theta) cos(Theta) 0; 0 0 1];
UzPrime = [1 0 0; 0 cos(Psi) -sin(Psi); 0 sin(Psi) cos(Psi)];
U = UzPrime * Ux * Uz;


end

