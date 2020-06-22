function [x_img, y_img, affine_mat] = latt2img_coord(transform_params, x_latt, y_latt)
% [x_img, y_img] = latticeToImg(x_latt, y_latt, P)
%
% P = [Theta1, Theta2, Phi1, Phi2, Lambda1, Lambda2, Offset1, Offset2, BorderPixX, BorderPixY]
%
% The numbering schemes correspond with lattice matrix coordinates. I.e.
% Theta1, Phi1, Lambda1 are associated with the Ylatt direction, since this
% would be the first coordinate in matrix indexing.
%
% Theta1 is the angle the 1-lattice axis makes with the horizontal. i.e.
% Theta1 = 90deg => 1-lattice direction is the same as y direction.
%
% Theta2 is defined in the same way. Everything is symmetric wrt to 1 and 2
%
% Lambda1 is the length in image space between lattice sites along lattice
% axis 1.
%
% Phi1 is the fractional number of lattice sites to offset the first index
% of the lattice coords by
%
% Offset1 is the integer number of lattice sites to offset the first index 
% of the lattice by.
%
% border_pix_x are the number of pixels at the border of the fluorescence
% image which are ignored by the reconstruction. This is an effective
% offset.
%
% Transformations
% [Ximg; Yimg] = [ [cos(Theta2) cos(Theta1)] ; [sin(Theta2) sin(Theta1)] ] * ...
%                [ [Lambda2 0]; [0 Lambda1] ] * 
%                ( [Xlatt; Ylatt] + [Offset2; Offset1] - [Phi2; Phi1] ) + ...
%                [border_pix_x; border_pix_y] 
%
% [Xlatt; Ylatt] = [ [1\lambda2 0]; [0 1/lambda1] ] * ...
%                  [ [-sin(Theta1) cos(Theta1)]; [sin(Theta2) -cos(Theta2)] ] * ...
%                  1/sin(Theta2-Theta1) * ( [Ximg; Yimg] - [border_pix_x; border_pix_y] ) - ...
%                  [Offset2; Offset1] + [Phi2; Phi1]
%
%TODO: Add clearer explanation of parameters...

theta1 = transform_params(1);
theta2 = transform_params(2);
phi1 = transform_params(3);
phi2 = transform_params(4);
lambda1 = transform_params(5);
lambda2 = transform_params(6);
offset1 = transform_params(7);
offset2 = transform_params(8);
border_pix_x = transform_params(9);
border_pix_y = transform_params(10);

% first translate
transl_vect = [-(phi2 - offset2); -(phi1 - offset1)];
[x_intermediate, y_intermediate, affine_xform1] = affine_transform([1 0; 0 1], transl_vect, x_latt, y_latt);
% do scaling transformation + rotation
xform_params = [lambda2 * sin(theta2) lambda1 * sin(theta1); lambda2 * cos(theta2) lambda1 * cos(theta1)];
[x_img, y_img, affine_xform2] = affine_transform(xform_params, [0; 0], x_intermediate, y_intermediate);
% finally, translate for border pix
border_transl_vect = [border_pix_x; border_pix_y];
[x_img, y_img, affine_xform3] = affine_transform([1 0; 0 1], border_transl_vect, x_img, y_img);

% kind of obnoxious that matlab uses affine transformation for row vectors
% instead of column vectors.
affine_mat = transpose(affine_xform1.T * affine_xform2.T * affine_xform3.T);

% x_img = (x_latt + offset2 - phi2) * lambda2 * sin(theta2) + ...
%         (y_latt + offset1 - phi1) * lambda1 * sin(theta1);
%    
% y_img = (x_latt + offset2 - phi2) * lambda2 * cos(theta2) + ....
%         (y_latt + offset1 - phi1) * lambda1 * cos(theta1);
    
% PSchauss implementation from reconstruction class.
% x = (indexCoord(1)+indexOffset(1)-g.phi1)*g.lambda(1)*cos(g.theta1) + ...
%     (indexCoord(2)+indexOffset(2)-g.phi2)*g.lambda(2)*cos(g.theta2);
% y = -(indexCoord(1)+indexOffset(1)-g.phi1)*g.lambda(1)*sin(g.theta1) - ...
%      (indexCoord(2)+indexOffset(2)-g.phi2)*g.lambda(2)*sin(g.theta2);
% xy = [ x -y ]; %
end