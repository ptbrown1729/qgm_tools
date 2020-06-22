function [img_d2x, unc_d2x, std_d2x, sem_d2x,...
          img_d2y, unc_d2y, std_d2y, sem_d2y] = ...
                                d2_avg_vectorfn(img_x, img_x_sem, img_y, img_y_sem)
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


img_d2x = 0.25 * d2_sum_x(img_x, img_y, 1);
% (img_x + rot90(img_y, 1) - rot90(img_x, 2) - rot90(img_y, 3));
        
img_d2y = 0.25 * d2_sum_y(img_x, img_y, 1);
% (img_y - rot90(img_x, 1) - rot90(img_y, 2) + rot90(img_x, 3));
        
        
% TODO: since some points are averaged with themselves, the uncertainty for
% those points should not be reduced.
% unc_x_sqr = unc_x.^2;
% unc_y_sqr = unc_y.^2;
% 
% unc_x_d2avg = 0.25 * sqrt(unc_x_sqr + rot90(unc_y_sqr, 1) + rot90(unc_x_sqr, 2) + rot90(unc_y_sqr, 3));
%         
% unc_y_d2avg = 0.25 * sqrt(unc_y_sqr + rot90(unc_x_sqr, 1) + rot90(unc_y_sqr, 2) + rot90(unc_x_sqr, 3));

std_d2x = sqrt( d2_sum_x( (img_x - img_d2x).^2, (img_y - img_d2y).^2, 0) / 3);
sem_d2x = std_d2x / sqrt(4);
unc_d2x = 1 ./ sqrt( d2_sum_x( 1 ./ img_x_sem.^2 , 1 ./ img_y_sem.^2, 0) ); 

std_d2y = sqrt( d2_sum_y( (img_x - img_d2x).^2, (img_y - img_d2y).^2, 0) / 3);
sem_d2y = std_d2y / sqrt(4);
unc_d2y = 1 ./ sqrt( d2_sum_y( 1 ./ img_x_sem.^2 , 1 ./ img_y_sem.^2, 0) ); 

end

function mat_d2x = d2_sum_x(mat_x, mat_y, signed)
    if signed == 1
        factor = -1;
    else
        factor = 1;
    end
    
    mat_d2x = mat_x ...
           + rot90(mat_y, 1) ...
           + factor *  rot90(mat_x, 2) ...
           + factor * rot90(mat_y, 3);
end

function mat_d2y = d2_sum_y(mat_x, mat_y, signed)
    if signed == 1
        factor = -1;
    else
        factor = 1;
    end
    
    mat_d2y = mat_y...
           + factor * rot90(mat_x, 1) ...
           + factor * rot90(mat_y, 2) ...
           + rot90(mat_x, 3);
end
