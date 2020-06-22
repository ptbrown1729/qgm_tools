function [V, PNames, ArgStr] = mask2d(P, X, Y, mask_fn)
%[V, PNames, ArgStr] = mask2D(P, X, Y)
%P = [Cx, Cy, A, Theta, Scale, Bg, Inversion]
PNames = {'Cx', 'Cy', 'A', 'Theta', 'Scale', 'Bg', 'Inversion'};
ArgStr = '';

%scaled,shifted,rotated,function



if isempty(P) || isempty(X) || isempty(Y)
    V = 0;
else
    [Xobj, Yobj, ~] = img2obj_coord(P, X, Y);
    V = mask_fn(Xobj, Yobj);
    V(isnan(V)) = 0;
end

end
