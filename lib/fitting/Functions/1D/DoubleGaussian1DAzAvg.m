function [V,PNames,ArgStr] = DoubleGaussian1DAzAvg(P,X)
%P = [SG1,AG1,SG2,AG2,Bg]
%1D Double Gaussian with centers at 0 and always > 0 amplitude
%The point of Fit1D and Fit2D is that you don't have to write
%combined functions like this. They can create the combined functions.
PNames = {'SG1','AG1','SG2','AG2','Bg'};
ArgStr = 'abs(AG1)*exp(-0.5*X.^2./(SG1^2)) + abs(AG2)*exp(-0.5*X.^2./(SG2^2))+ Bg';

if isempty(P)||isempty(X)
    V = 0;
else
    SG1 = P(1); AG1 = P(2);
    SG2 = P(3); AG2 = P(4); Bg = P(5);
    V = abs(AG1)*exp(-0.5*X.^2./(SG1^2)) + abs(AG2)*exp(-0.5*X.^2./(SG2^2))+ Bg;
end
end
