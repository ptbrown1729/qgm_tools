function [V,PNames,ArgStr] = exp1D(P,X)
%[V,PNames,ArgStr] = exp1D(P,X)
%P = [Cx,A,Tau,Bg]
%A*exp(-(X-Cx)./Tau)+Bg;
PNames = {'Cx','A','Tau','Bg'};
ArgStr = 'A*exp((X-Cx)./Tau)+Bg';

if isempty(P)||isempty(X)
    V = 0;
else
    Cx = P(1); A = P(2); Tau = P(3); Bg = P(4);
    if Tau~=0
        V = A*exp(-(X-Cx)./Tau)+Bg;
    else
        V = zeros(size(X));
    end
end
end
