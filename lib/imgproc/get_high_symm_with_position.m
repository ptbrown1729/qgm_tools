function [array_parts, x_parts, y_parts] = get_high_symm_with_position(array_gxmg, x_gxmg, y_gxmg, linear_index)


% GX cut
indices_gx = linear_index <= 1 & linear_index >= 0;
array_gx = array_gxmg(indices_gx);
x_gx = x_gxmg(indices_gx);
y_gx = y_gxmg(indices_gx);

% XM cut
indices_xm = linear_index >= 1 & linear_index <= 2;
array_xm = array_gxmg(indices_xm);
x_xm = x_gxmg(indices_xm);
y_xm = y_gxmg(indices_xm);

% MG cut
indices_mg = linear_index >= 2 & linear_index <= 3;
array_mg = array_gxmg(indices_mg);
x_mg = x_gxmg(indices_mg);
y_mg = y_gxmg(indices_mg);

%
array_parts = {array_gx, array_xm, array_mg};
x_parts = {x_gx, x_xm, x_mg};
y_parts = {y_gx, y_xm, y_mg};

end