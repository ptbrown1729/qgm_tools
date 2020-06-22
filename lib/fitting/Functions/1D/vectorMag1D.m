function [V,PNames,ArgStr] = vectorMag1D(P,X)
%[V] = vectorMag1D(P,X,PNamesBool)
%V = A*sqrt(B^2+(x-xc)^2)+bg
%P = [A,B,xc,bg]
%1D Line
PNames = {'Amp','By','xc','bg'};
ArgStr = 'Amp*sqrt(By^2+(X-xc).^2)+bg';

if isempty(P)||isempty(X)
    V = 0;
else
    Amp = P(1);
    By = P(2);
    xc = P(3);
    bg = P(4);
    V = Amp * sqrt(By^2 + (X - xc) .^2) + bg;
end
end


