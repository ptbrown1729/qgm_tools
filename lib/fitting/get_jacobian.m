function [jacobian_mat, f_init_pt] = get_jacobian(fn, init_point, steps)

jacobian_mat = zeros(1, length(init_point));

f_init_pt = fn(init_point);
for jj = 1:length(init_point)
    new_pt = init_point;
    new_pt(jj) = new_pt(jj) + steps(jj);
	jacobian_mat(jj) = (fn(new_pt) - f_init_pt) / (steps(jj));
end

end