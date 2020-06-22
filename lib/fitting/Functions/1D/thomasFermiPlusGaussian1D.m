function [V,PNames,ArgStr] = thomasFermiPlusGaussian1D(P,X)
%function [V] = thomasFermiPlusGaussian1D(P,X,PNamesBool)
%P = [CTF,RTF,ATF,CG,SG,AG,Bg]
%1D Thomas-Fermi plus Gaussian
%The point of Fit1D and Fit2D is that you don't have to write
%combined functions like this. fit1D can combine the functions for you.
PNames = {'CTF','RTF','ATF','CG','SG','AG','Bg'};
ArgStr = '';

if isempty(P)||isempty(X)
    V = 0;
else
    CTF = P(1); RTF = P(2); ATF = P(3);
    CG = P(4); SG = P(5); AG = P(6); Bg = P(7);
    Vtemp = ATF*(1-(X-CTF).^2/RTF^2);
    V =Vtemp.*(Vtemp>0) + AG*exp(-0.5*(X-CG).^2./(SG^2))+ Bg;
end
end
