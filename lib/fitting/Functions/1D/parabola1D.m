function [V, PNames, ArgStr] = parabola1D(P, X)
%[V] = parabola1D(P,X,PNamesBool)
%P = [P1,P2,P3]
%P1*X.^2+P2*X+P3
PNames = {'P1','P2','P3'};
ArgStr = 'P1*X.^2+P2*X+P3';

if isempty(P) || isempty(X)
    V = 0;
else
    P1 = P(1);
    P2 = P(2);
    P3 = P(3);
    V = P1 * X .^ 2 + P2 * X + P3;
end
end


