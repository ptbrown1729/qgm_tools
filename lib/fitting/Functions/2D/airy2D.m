function [V,PNames,ArgStr] = airy2D(P,X,Y)
%[V,PNames,ArgStr] = airy2D(P,X,Y)
%P = [Cx,Cy,A,FirstZero,Bg]

PNames = {'Cx','Cy','A','FirstZero','Bg'};
ArgStr = '';

if isempty(P)||isempty(X)
    V = 0;
else
    Cx = P(1); 
    Cy = P(2); 
    A = P(3); 
    Rzero = P(4); 
    Bg = P(5);
    
    ScaleFact = 3.8317 / Rzero;
    R = sqrt((X - Cx).^2 + (Y - Cy).^2);
    V = A * (2 * besselj(1, R * ScaleFact) ./ (R * ScaleFact)).^2 + Bg;
    V(R == 0) = A + Bg;
end

end
