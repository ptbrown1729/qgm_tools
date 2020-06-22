function [x_img, y_img, affine_mat] = obj2img(transform_params, x_obj, y_obj)
% Perform a specific type of affine transformation on coordinates. We can 
% regard this a transformation between 'object' and 'image' space, but it
% is generally a coordinate transformation where we allow shifting of the
% origin, rotation, and dilation.
%
% transform_params = [cx_img, cy_img, theta, scale] are the parameters of the transformation
% The transformation parameters are defined so that img2obj is the inverse
% transformation for obj2img if the same parameters are used.%
%
% When transforming from object space to image space:
% Theta convention. Can regard this as rotating vectors CW, or the
% coordinate system CCW.
%
% If Scale > 1, vectors are expanded, or the coordinate system is
% contracted
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

transform_mat = scale * [cos(theta) sin(theta); -sin(theta) cos(theta)];
[x_img, y_img, affine_xform1] = affine_transform(transform_mat, [cx_img; cy_img], x_obj, y_obj);

% x_img = x_obj * cos(theta) * scale + y_obj * sin(theta) * scale + cx_img;
% y_img = -x_obj * sin(theta) * scale + y_obj * cos(theta) * scale + cy_img;

% kind of obnoxious that matlab uses affine transformation for row vectors
% instead of column vectors.
affine_mat = transpose(affine_xform1.T);

end