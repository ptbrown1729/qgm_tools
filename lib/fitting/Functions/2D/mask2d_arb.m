function [V, PNames, ArgStr] = mask2d_arb(P, X, Y, mask, x_mask, y_mask)
%[V, PNames, ArgStr] = mask2D(P, X, Y)
%P = [Cx, Cy, Theta, Scale]
PNames = {'Cx', 'Cy', 'Theta', 'Scale'};
ArgStr = '';

if isempty(P) || isempty(X)
    V = 0;
else
    Cx = P(1); 
    Cy = P(2); 
    Theta = P(3); 
    Scale = P(4); 
    
    % X, Y are in image space and mask is in the object space
    
    if ~exist('x_mask', 'var') || ~exist('y_mask', 'var')
        [x_mask, y_mask] = meshgrid(1:size(mask, 2), 1:size(mask, 1));
    end
    %scaled,shifted,rotated,function
    %[cx_img, cy_img, theta, scale]
    xform_params = [Cx, Cy, Theta, Scale]; 
    % Given image space point [x_i, y_i], what are the corresponding object
    % space coordinates [x_o, y_o]
    [x_obj, y_obj] = img2obj_coord(xform_params, X, Y);

    % assuming that object coordinates are 1, ... size of mask
    % we want to evaluate our mask at [x_o, y_o], M(y_o, x_o). Since we've
    % divided up our image space into pixels, we take the closest pixel to
    % (y_o, x_o).

 
%     x_indices = round(x_obj);
%     y_indices = round(y_obj);
% 
%     index_mask = ones(size(x_indices));
%     index_mask(x_obj > max(max(x_mask))) = 0;
%     index_mask(x_obj < min(min(x_mask))) = 0;
%     index_mask(y_obj > max(max(y_mask))) = 0;
%     index_mask(y_obj < min(min(y_mask))) = 0;
%     x_indices(index_mask == 0) = 1;
%     y_indices(index_mask == 0) = 1;
% 
%     mask_img = zeros(size(X));
%     mask_img(:, :) = mask(sub2ind(size(mask), y_indices, x_indices));
%     mask_img(index_mask == 0) = 0;
%     
%     V = mask_img;
    
    % alternative approach: create an interpolating function
    mask_interp = @(xo, yo) interp2(x_mask, y_mask, mask, xo, yo);
    mask_img = mask_interp(x_obj, y_obj);
    V = mask_img;
    V(isnan(V)) = 0;
end

end
