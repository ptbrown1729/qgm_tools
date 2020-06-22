function [V,PNames,ArgStr] = dampedSinusoid1D(P,X)
%[V,PNames,ArgStr] = dampedSinusoid1D(P,X,PNamesBool)
%P = [Cx,F,A,Phase,Bg,Tau]
PNames = {'Cx','F','A','Phase','Bg','Tau'};
ArgStr = 'A*sin(2*pi*F*(X-Cx)+Phase).*exp(-(X-Cx)/Tau)+Bg';

if isempty(P)||isempty(X)
    V = zeros(size(X));
else
    Cx = P(1); F = P(2); A = P(3); Phase = P(4); Bg = P(5); Tau = P(6);
    if Tau == 0
        V = zeros(size(X));
    else
        V = A*sin(2*pi*F.*(X-Cx)+Phase).*exp(-(X-Cx)./Tau)+Bg;
    end
end

end



