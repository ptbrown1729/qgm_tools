function [V,PNames,ArgStr] = line1D(P,X)
%[V] = line1D(P,X,PNamesBool)
%P = [Slope,Intercept]
%1D Line
PNames = {'slope','intercept'};
ArgStr = 'slope*X+intercept';

if isempty(P)||isempty(X)
    V = 0;
else
    slope = P(1); intercept = P(2);
    V = slope*X+intercept;
end

end
