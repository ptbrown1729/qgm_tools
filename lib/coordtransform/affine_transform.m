function [xs_transformed, ys_transformed, affine_tform] = affine_transform(linear_transform, transl_vector, xs, ys)

% actually, matlab has affine2d objects, which might be a better option
% than writing this from scratch. Except affine2d uses the matrix that
% would be applied to a row vector instead of a column vector as input
% argument. That is morally objectionable, so maybe I will just wrap their
% method with this...

ndim = size(linear_transform, 1);
affine_mat = zeros(ndim + 1, ndim + 1);

affine_mat(1:ndim, 1:ndim) = linear_transform;
affine_mat(end, 1:ndim) = 0;
affine_mat(1:ndim, end) = transl_vector(:);
affine_mat(end, end) = 1;

% coords = zeros(3, numel(xs));
% coords(1, :) = xs(:);
% coords(2, :) = ys(:);
% coords(3, :) = 1;
% 
% coords_transformed = affine_mat * coords;
% xs_transformed = coords_transformed(1, :);
% ys_transformed = coords_transformed(2, :);
% 
% xs_transformed = reshape(xs_transformed, size(xs));
% ys_transformed = reshape(ys_transformed, size(ys));

affine_tform = affine2d(transpose(affine_mat));
[xs_transformed, ys_transformed] = transformPointsForward(affine_tform, xs, ys);

end