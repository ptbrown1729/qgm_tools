function [V,PNames,ArgStr] = gaussianBeamWaist1D(P,X)
%[V,PNames,ArgStr] = gaussianBeamWaist1D(P,X)
%P = [Cx,W0,Lambda]
PNames = {'Cx','W0','Lambda'};
ArgStr = 'W0.*sqrt(1+((Cx-X)./(W0.^2*pi./Lambda.^2)).^2)';

if isempty(P)||isempty(X)
    V = zeros(size(X));
else
    Cx = P(1); W0 = P(2); Lambda = P(3);
    V = W0.*sqrt(1+((X-Cx)./(W0.^2*pi./Lambda)).^2);
end

end



