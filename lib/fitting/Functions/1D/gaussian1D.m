function [V,PNames,ArgStr] = gaussian1D(P,X)
%[V,PNames,ArgStr] = gaussian1D(P,X,PNamesBool)
%P = [CG,SG,AG,Bg]
%NOTE: None of the arguments can be single ascii letters.
PNames = {'CG','SG','AG','Bg'};
ArgStr = 'AG*exp(-0.5*(x-CG).^2/SG^2) + Bg';

if isempty(P)||isempty(X)
    V = 0;
else
    CG = P(1);SG = P(2);AG = P(3);Bg = P(4);
    V = AG*exp(-0.5*(X-CG).^2./(SG^2))+Bg;
end

end

