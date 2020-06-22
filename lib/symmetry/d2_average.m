function [img_d2_avg, unc_d2, img_d2_std, img_d2_sem] = d2_average(img, img_sem)
% [img_d2_avg, unc_d2, img_d2_std, img_d2_sem] = d2_average(img, img_sem)
% d2_average Produces an average of a matrix assuming that matrix possess
% rectangular symmetry. For example, the points at all four corners will be
% averaged together. We refer to this as 'd2 averaging' because the
% symmetry group of the rectangle is the second dihedral group. It has 4
% symmetry operations: the identity, 180 degree rotation, reflection, and
% 180 degree rotation combined with reflection
%
%   arguments:
%   -----------------
%   img: the 2D matrix to be average
%
%   img_sem: the standard errors of the mean associated with each point in
%   the image
%
%   outputs:
%   -----------------
%   img_d2_avg: the matrix average being D2 averaged. It is the same size
%   as the original matrix and therefore contains redundant information. A
%   single quadrant of this matrix contains the complete information.
%   
%   unc_d2: the weighted uncertainties associated with each pixel of
%   img_d2_avg. These are computed according to unc_I = 1 / sqrt( sum( 1 / sigma_i^2 ) )
%   where I refers to a pixel in the final image, and the i run over the
%   four pixels which are averaged together
%
%   img_d2_std: The standard deviation of the four pixels which are being
%   averaged. This includes 'Bessels correction', i.e. we divide by 1/N-1 
%   instead of 1/N
%
%   img_d2_sem: The standard error of the mean of the four pixels which are
%   bing averaged. This is equal to img_d2_std / sqrt(4).

if ~exist('img_sem', 'var') || isempty(img_sem)
    img_sem = zeros(size(img));
end

img_d2_avg = 0.25 * d2_sum(img);

% TODO: since some points are averaged with themselves, the uncertainty for
% those points should not be reduced.
% standard deviation and standard error of the mean
npix = 4;
img_d2_std = sqrt( d2_sum( (img_d2_avg - img).^2) / (npix - 1) );
img_d2_sem = img_d2_std / sqrt(npix);

% weighted uncertainty
% ignore points with zero uncertainty for calculating uncertainty
unc_weight = 1 ./ img_sem .^2;
unc_weight( isinf(unc_weight) ) = 0;

unc_d2 =  1 ./ sqrt( d2_sum( unc_weight ) );

% handle case where none of the points in the ROI have finite uncertainty.
unc_d2( isinf(unc_d2 )) = 0;

end

function d2_mat_sum = d2_sum(mat)
    d2_mat_sum = (mat + rot90(mat, 2) + rot90(transpose(mat), 1) + rot90(transpose(mat), 3));
end
