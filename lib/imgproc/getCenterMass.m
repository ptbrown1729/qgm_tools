function [Cx,Cy,Sx,Sy,Weight] = getCenterMass(M)
%[Cx,Cy,Sx,Sy,Weight] = getCenterMass(M)
%Accepts a stack of images.
% TODO: deprecate in favor of get_moment

[SizeY,SizeX,Nimgs] = size(M);
Cx = zeros(Nimgs,1);
Cy = zeros(Nimgs,1);
Sx = zeros(Nimgs,1);
Sy = zeros(Nimgs,1);
Weight = zeros(Nimgs,1);
for ii = 1:Nimgs
    try
    CurrentM = M(:,:,ii);
    Vx = sum(CurrentM,1);
    Weight(ii) = sum(Vx);
    Vy = sum(CurrentM,2)';

    %Vx = Vx.*(Vx>0);
    %Vy = Vy.*(Vy>0);

    X = 1:SizeX;
    Y = 1:SizeY;

    Cx(ii) = sum(Vx.*X)/sum(Vx);
    Cy(ii) = sum(Vy.*Y)/sum(Vy);

    Sx(ii) = sqrt(abs(sum(Vx.*(abs(X-Cx(ii)).^2))/sum(Vx)));
    Sy(ii) = sqrt(abs(sum(Vy.*(abs(Y-Cy(ii)).^2))/sum(Vy)));
    catch
        Cx(ii) = NaN;
        Cy(ii) = NaN;
        Sx(ii) = NaN;
        Sy(ii) = NaN;
    end
end

end
