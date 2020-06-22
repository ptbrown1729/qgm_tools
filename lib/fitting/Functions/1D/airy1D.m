function [V,PNames,ArgStr] = airy1D(P,X)
%[V,PNames,ArgStr] = airy1D(P,X)
%P = [Cx,A,FirstZero,Bg]
PNames = {'Cx','A','FirstZero','Bg'};
ArgStr = 'A*(2*besselj(1,(X-Cx)*ScaleFact)./((X-Cx)*ScaleFact)).^2+Bg';

if isempty(P)||isempty(X)
    V = 0;
else
    Cx = P(1); A = P(2); Xzero = P(3); Bg = P(4);
    ScaleFact = 3.8317/Xzero;
    V = A*(2*besselj(1,(X-Cx)*ScaleFact)./((X-Cx)*ScaleFact)).^2+Bg;
    V(X==Cx) = A + Bg;
    

end

end
