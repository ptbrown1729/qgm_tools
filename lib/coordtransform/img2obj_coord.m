function [x_obj, y_obj, affine_mat] = img2obj_coord(transform_params, x_img, y_img)
% Perform a specific type of affine transformation on coordinates. We can 
% regard this a transformation between 'object' and 'image' space, but it
% is generally a coordinate transformation where we allow shifting of the
% origin, rotation, and dilation.
%
% transform_params = [cx_img, cy_img, theta, scale] are the parameters of the transformation
% The transformation parameters are defined so that obj2img is the inverse
% transformation for img2obj if the same parameters are used.
% 
%
% When transforming from image space to object space:
% Theta convention. Can regard this as rotating vectors CCW, or the
% coordinate system CW.
%
% If Scale > 1, vectors are contracted, or the coordinate system is
% expanded.
%
% The transformation and inverse transformation are given by...
% [Xobj; Yobj] = [[cos(Theta) -sin(Theta)]; [sin(Theta) cos(Theta)]] * ([Ximg; Yimg]-[Cx_img; Cy_img]) * 1/Scale
% [Ximg; Yimg] = Scale * [[cos(Theta) sin(Theta)]; [-sin(Theta) cos(Theta)]] * [Xobj; Yobj] + [Cx_img; Cy_img]
%              = Scale * [[cos(Theta) sin(Theta)]; [-sin(Theta) cos(Theta)]] * [Xobj; Yobj] + [Cx_img; Cy_img]
%
% Given a function f(xi, yi) defined on image space, and a transformation
% xi = T1(xo,yo), yi =T2(xo,yo), then we define the same function on the
% object space according to fo(xo, yo) = f(T1(xo,yo), T2(xo,yo)). It is
% more natural to write the inverse transformation which converts
% coordinates from image space to object space, so it is more standard to
% see f(inv(T)(xi, yi))

cx_img = transform_params(1);
cy_img = transform_params(2);
theta = transform_params(3);
scale = transform_params(4);

% first translate
[x_intermediate, y_intermediate, affine_xform1] = affine_transform([1 0; 0 1], [-cx_img; -cy_img], x_img, y_img);
% then scale/rotate
transform_mat = 1 / scale * [cos(theta) -sin(theta); sin(theta) cos(theta)];
[x_obj, y_obj, affine_xform2] = affine_transform(transform_mat, [0; 0], x_intermediate, y_intermediate);

% x_obj = (x_img-cx_img) * cos(theta) / scale - (y_img - cy_img) * sin(theta) / scale;
% y_obj = (x_img-cx_img) * sin(theta) / scale + (y_img - cy_img) * cos(theta) / scale;

% kind of obnoxious that matlab uses affine transformation for row vectors
% instead of column vectors.
affine_mat = transpose(affine_xform1.T * affine_xform2.T);

end
