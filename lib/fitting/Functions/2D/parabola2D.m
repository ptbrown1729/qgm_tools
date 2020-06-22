function [val, param_names, arg_str] = parabola2D(P, X, Y)
% Rotated 2D Parabola
%
% [V] = [V, PNames, ArgStr] = parabola2D(P,X,Y)
% P = [cx, cy, ax, ay, theta, bg]
% Theta is the counter-clockwise angle between the axis with width SxG, if
% we regard the y axis as pointing downards (i.e. as matlab would display a
% matrix).

param_names = {'cx', 'cy', 'ax', 'ay', 'theta', 'Bg'};
arg_str = '';

if isempty(P) || isempty(X) || isempty(Y)
    val = 0;
else
    cx = P(1); 
    cy = P(2); 
    ax = P(3); 
    ay = P(4); 
    theta = P(5); 
    bg = P(6);
    
    val = bg + ax * ( (X - cx) * cos(theta) - (Y - cy) * sin(theta)).^2 ...
           + ay * ( (Y - cy) * cos(theta) + (X - cx) * sin(theta)).^2;
end
end

