function [img_d4x, unc_d4x, std_d4x, sem_d4x,...
          img_d4y, und_d4y, std_d4y, sem_d4y] = ...
                                d4_avg_vectorfn(img_x, img_x_sem, img_y, img_y_sem)
%D4_AVERAGE Summary of this function goes here
%   Detailed explanation goes here
%   Given a vector function V(x,y) = [f(x,y), g(x,y)], the transformation law is 
%   Vt(x,y) = T V( T^{-1}(x,y) )
%   
%   Need to think if I care about this here?ar
%
%   For example, consider the 90 degree rotation matrix 
%   T^-1 = [0 1; -1 0]
%   T    = [0 -1; 1 0]
%   Vt(x,y) = [ -g(y, -x), f(y, -x) ]
%
%   TODO: finish ... do I really need these minus signs?
%
if ~isequal( size(img_x), size(img_y))
    error();
end

if ~exist('unc_x', 'var') || isempty(img_x_sem)
    img_x_sem = zeros(size(img_x));
end

if ~exist('unc_y', 'var') || isempty(img_y_sem)
    img_y_sem = zeros(size(img_x));
end


img_d4x = 0.125 * d4_sum_x(img_x, img_y, 1);
% 0.125 * (img_x + rot90(img_y, 1) - rot90(img_x, 2) - rot90(img_y, 3) ...
%             + transpose(img_y) + rot90(transpose(img_x), 1) ...
%             - rot90(transpose(img_y), 2) - rot90(transpose(img_x), 3));
        
img_d4y = 0.125 * d4_sum_y(img_x, img_y, 1);
% 0.125 * (img_y - rot90(img_x, 1) - rot90(img_y, 2) + rot90(img_x, 3) ...
%             + transpose(img_x) - rot90(transpose(img_y), 1) ...
%             - rot90(transpose(img_x), 2) + rot90(transpose(img_y), 3) );
        
        
% TODO: since some points are averaged with themselves, the uncertainty for
% those points should not be reduced.

% weighted average 

std_d4x = sqrt( d4_sum_x( (img_x - img_d4x).^2, (img_y - img_d4y).^2, 0) / 7);
sem_d4x = std_d4x / sqrt(8);
unc_d4x = 1 ./ sqrt( d4_sum_x( 1 ./ img_x_sem.^2 , 1 ./ img_y_sem.^2, 0) ); 
% 0.125 * sqrt(unc_x_sqr + rot90(unc_y_sqr, 1) + rot90(unc_x_sqr, 2) + rot90(unc_y_sqr, 3) + ...
%             transpose(unc_y_sqr) + rot90(transpose(unc_x_sqr), 1) + ...
%             rot90(transpose(unc_y_sqr), 2) + rot90(transpose(unc_x_sqr), 3));
        
% unc_y_d4avg = 
% 0.125 * sqrt(unc_y_sqr + rot90(unc_x_sqr, 1) + rot90(unc_y_sqr, 2) + rot90(unc_x_sqr, 3) + ...
%             transpose(unc_x_sqr) + rot90(transpose(unc_y_sqr), 1) + ...
%             rot90(transpose(unc_x_sqr), 2) + rot90(transpose(unc_y_sqr), 3));
std_d4y = sqrt( d4_sum_y( (img_x - img_d4x).^2, (img_y - img_d4y).^2, 0) / 7);
sem_d4y = std_d4y / sqrt(8);
und_d4y = 1 ./ sqrt( d4_sum_y( 1 ./ img_x_sem.^2 , 1 ./ img_y_sem.^2, 0) ); 

end

function mat_d4x = d4_sum_x(mat_x, mat_y, signed)
    if signed == 1
        factor = -1;
    else
        factor = 1;
    end
    
    mat_d4x = mat_x ...
           + rot90(mat_y, 1) ...
           + factor *  rot90(mat_x, 2) ...
           + factor * rot90(mat_y, 3) ...
           + transpose(mat_y) ...
           + rot90(transpose(mat_x), 1) ...
           + factor * rot90(transpose(mat_y), 2) ...
           + factor * rot90(transpose(mat_x), 3);
end

function mat_d4y = d4_sum_y(mat_x, mat_y, signed)
    if signed == 1
        factor = -1;
    else
        factor = 1;
    end
    
    mat_d4y = mat_y...
           + factor * rot90(mat_x, 1) ...
           + factor * rot90(mat_y, 2) ...
           + rot90(mat_x, 3) ...
           + transpose(mat_x)...
           + factor * rot90(transpose(mat_y), 1) ...
           + factor * rot90(transpose(mat_x), 2) ...
           + rot90(transpose(mat_y), 3);
end
