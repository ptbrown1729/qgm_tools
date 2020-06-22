function [x_latt, y_latt, affine_mat] = img2latt_coord(transform_params, x_img, y_img)
% [x_latt, y_latt] = imgToLattice(Ximg,Yimg,P)
%
% P = [Theta1, Theta2, Phi1, Phi2, Lambda1, Lambda2, Offset1, Offset2]
%
% The numbering schemes correspond with lattice matrix coordinates. I.e.
% Theta1, Phi1, Lambda1 are associated with the Ylatt direction, since this
% would be the first coordinate in matrix indexing.
%
% Theta1 is the angle the 1-lattice axis makes with the horizontal. i.e.
% Theta1 = 90deg => 1-lattice direction is the same as y direction.
%
% Theta2 is defined in the same way. Everything is symmetric wrt to Theta1
% and Theta2
%
% Lambda1 is the length in image space between lattice sites along lattice
% axis 1.
%
% Offset1 is the number of lattice sites to offset the first index of the
% lattice by.
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

%Ximg and Yimg should probably be restricted to integers...
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

% first account for border pixels
border_transl_vect = [-border_pix_x; -border_pix_y];
[x_img, y_img, affine_xform0] = affine_transform([1 0; 0 1], border_transl_vect, x_img, y_img);
% first rotate to get the lattice axes
% theta1 = 0, the new axis lies along the -x direction. theta1 is positive
% going CW from -x
% theta2 = 0, the new axis lies along the +x direction. theta2 is positive
% going CW from x.
xform_params = [-cos(theta1) sin(theta1); cos(theta2) -sin(theta2)];
[x_intermediate, y_intermediate, affine_xform1] = affine_transform(xform_params, [0; 0], x_img, y_img);
% do scaling transformation + offset
xform_params = 1 / sin(theta1 - theta2) * [1 / lambda2 0; 0 1 /lambda1];
transl_vect = [phi2 - offset2; phi1 - offset1];
[x_latt, y_latt, affine_xform2] = affine_transform(xform_params, transl_vect, x_intermediate, y_intermediate);

% kind of obnoxious that matlab uses affine transformation for row vectors
% instead of column vectors.
affine_mat = transpose(affine_xform0.T * affine_xform1.T * affine_xform2.T);

% x_latt = ( -x_img * cos(theta1) + y_img * sin(theta1) ) * 1/lambda2 * 1 / sin(theta1 - theta2) + phi2 - offset2;
% y_latt = ( x_img * cos(theta2) - y_img * sin(theta2) ) *1/lambda1 * 1 / sin(theta1 - theta2) + phi1 - offset1;

% x(2)=-x(2);% added "-" to fit old trafo
% t = g.transformFactor;
% n = (-x(1)*sin(g.theta2)-x(2)*cos(g.theta2))*t(1) + g.phi1;
% m = ( x(1)*sin(g.theta1)+x(2)*cos(g.theta1))*t(2) + g.phi2;
%
% nm = [round(n)-indexOffset(1) round(m)-indexOffset(2)];


end