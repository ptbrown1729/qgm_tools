function [V,PNames,ArgStr] = lorentzian1D(P,X)
%function [V] = lorentzian1D(P,X,PNamesBool)
%P = [Cx,HWHM,A,Bg]
%Function Name
PNames = {'Cx','HWHM','Amp','Bg'};
ArgStr = 'Amp./(1+((X-Cx)./HWHM).^2)+Bg';

if isempty(P)||isempty(X)
    V = 0;
else
    Cx = P(1); HWHM = P(2); Amp = P(3); Bg = P(4);
    V = Amp./(1+((X-Cx)./HWHM).^2)+Bg;
end
end



