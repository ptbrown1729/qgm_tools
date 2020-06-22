function [V,PNames,ArgStr] = TFGaussian1DAzAvg(P,X)
%P = [RTF,ATF,SG,AG,Bg]
%1D Thomas-Fermi plus Gaussian with centers at 0 and always > 0 amplitude
%The point of Fit1D and Fit2D is that you don't have to write
%combined functions like this. They can create the combined functions.
PNames = {'RTF','ATF','SG','AG','Bg'};
ArgStr = '';

if isempty(P)||isempty(X)
    V = 0;
else
    RTF = P(1); ATF = P(2);
    SG = P(3); AG = P(4); Bg = P(5);
    Vtemp = abs(ATF)*(1-X.^2./RTF^2);
    V =Vtemp.*(Vtemp>0) + abs(AG)*exp(-0.5*X.^2./(SG^2))+ Bg;
end
end
