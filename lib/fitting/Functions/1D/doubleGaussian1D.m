function [V,PNames,ArgStr] = doubleGaussian1D(P,X)
%P = [CG1,SG1,AG1,CG2,SG2,AG2,Bg]
%sum of two 1D Gaussians
%The point of fit1D is that it can combine functions for you...
PNames = {'CG1','SG1','AG1','CG2','SG2','AG2','Bg'};
ArgStr = 'AG1*exp(-0.5*(X-CG1).^2./(SG1^2))+AG2*exp(-0.5*(X-CG2).^2./(SG2^2))+Bg';

if isempty(P)||isempty(X)
    V = 0;
else
    CG1 = P(1);SG1 = P(2);AG1 = P(3);CG2 = P(4);SG2 = P(5);AG2 = P(6);Bg = P(7);
    V = AG1*exp(-0.5*(X-CG1).^2./(SG1^2))+AG2*exp(-0.5*(X-CG2).^2./(SG2^2))+Bg;
end
end
