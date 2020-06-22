function [V,PNames,ArgStr] = thomasFermi1D(P,X)
%function [V] = thomasFermi1D(P,X,PNamesBool)
%P = [CTF,RTF,ATF,Bg]
%1D Thomas-Fermi Function.
PNames = {'CTF','RTF','ATF','Bg'};
ArgStr = '';

if isempty(P)||isempty(X)
    V = 0;
else
    CTF = P(1); RTF = P(2); ATF = P(3); Bg = P(4);
    Vtemp = ATF*(1-(X-CTF).^2/RTF^2);
    V =Vtemp.*(Vtemp>0) + Bg;
end
end

